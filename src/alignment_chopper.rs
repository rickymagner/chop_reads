use std::cmp::min;
use rust_htslib::bam::{Record};
use rust_htslib::bam::record::{CigarString, Cigar, Aux};

pub struct AlignmentChopper {
    chunk_size: usize,
    include_edges: bool,
    skip_softclips: bool,
    read_group: Option<String>,
}

impl AlignmentChopper {
    pub fn new(chunk_size: usize, include_edges: bool, skip_softclips: bool, read_group: Option<String>) -> Self {
        Self {
            chunk_size,
            include_edges,
            skip_softclips,
            read_group
        }
    }

    pub fn chop_read(&self, rec: &Record) -> Vec<Record> {
        let num_pieces = if self.include_edges {
            (rec.seq_len() + self.chunk_size - 1) / self.chunk_size  // Rounds up
        } else {
            rec.seq_len() / self.chunk_size
        };

        // Take rev direction cigar string so can use Vec & pop as stack later
        // Note: Need to keep track of leading/trailing softclips to offset the start/end of POS for chunks
        let mut rev_cigar = if self.skip_softclips {
            let mut cigar_iter = rec.cigar().iter();
            if rec.cigar().leading_softclips() > 0 {
                cigar_iter.next();
            }
            let cigar_vec: Vec<Cigar> = if rec.cigar().trailing_softclips() > 0 {
                cigar_iter.rev().next().collect()
            } else {
                cigar_iter.collect()
            };
            CigarString(cigar_vec)
        } else {
            CigarString(rec.cigar().take().0.into_iter().rev().collect())
        };

        let mut rec_pieces = Vec::with_capacity(num_pieces);
        let mut ref_pos_offset = 0;
        for i in 0..num_pieces {
            let mut new_rec = Record::default();
            let chunk_begin = self.chunk_size * i;
            let chunk_end = min(self.chunk_size * (i+1), rec.seq_len());

            // Generate the new fields for the current chop slice
            let (new_cigar, num_ref_bases_consumed) = self.chop_cigar(&mut rev_cigar);
            let new_pos = rec.pos() + ref_pos_offset as i64;  // Use previous offset then update for next read
            ref_pos_offset += num_ref_bases_consumed;

            // Get seq and qual slices
            let new_seq = &rec.seq().as_bytes()[chunk_begin..chunk_end];
            let new_qual = &rec.qual()[chunk_begin..chunk_end];

            // Update name for chunk
            let new_qname = &[rec.qname(), b"-", i.to_string().as_bytes()].concat();

            // These are changed based on the particular slice
            new_rec.set(new_qname, Some(&new_cigar), new_seq, new_qual);
            new_rec.set_pos(new_pos);

            // Following are unchanged
            new_rec.set_flags(rec.flags());
            new_rec.set_tid(rec.tid());
            new_rec.set_mapq(rec.mapq());
            new_rec.set_mtid(rec.mtid());
            new_rec.set_mpos(rec.mpos());
            new_rec.set_insert_size(rec.insert_size());

            // All aux data other than RG is lost
            if let Some(rg) = &self.read_group {
                if let Ok(_a) = new_rec.aux(b"RG") {
                    new_rec.remove_aux(b"RG").expect(&format!("Could not remove RG from: {} - {}", &new_rec.tid(), &new_rec.pos()));
                }
                new_rec.push_aux(b"RG", Aux::String(rg)).expect(&format!("Unable to push RG string at: {} - {}", &new_rec.tid(), &new_rec.pos()));
            }

            rec_pieces.push(new_rec);
        }

        rec_pieces
    }

    fn consumes_query(c: &Cigar) -> bool {
        match c {
            Cigar::Match(_) => true,
            Cigar::Ins(_) => true,
            Cigar::Del(_) => false,
            Cigar::RefSkip(_) => false,
            Cigar::SoftClip(_) => true,
            Cigar::HardClip(_) => false,
            Cigar::Pad(_) => false,
            Cigar::Equal(_) => true,
            Cigar::Diff(_) => true
        }
    }

    fn get_num_ref_bases_consumed(c: &Cigar) -> u32 {
        match c {
            Cigar::Match(i) => *i,
            Cigar::Ins(_) => 0,
            Cigar::Del(i) => *i,
            Cigar::RefSkip(i) => *i,
            Cigar::SoftClip(_) => 0,
            Cigar::HardClip(_) => 0,
            Cigar::Pad(_) => 0,
            Cigar::Equal(i) => *i,
            Cigar::Diff(i) => *i
        }
    }

    fn split_cigar(c: &Cigar, pos: u32) -> (Cigar, Cigar) {
        match c {
            Cigar::Match(i) => (Cigar::Match(pos), Cigar::Match(*i-pos)),
            Cigar::Ins(i) => (Cigar::Ins(pos), Cigar::Ins(*i-pos)),
            Cigar::Del(i) => (Cigar::Del(pos), Cigar::Del(*i-pos)),
            Cigar::RefSkip(i) => (Cigar::RefSkip(pos), Cigar::RefSkip(*i-pos)),
            Cigar::SoftClip(i) => (Cigar::SoftClip(pos), Cigar::SoftClip(*i-pos)),
            Cigar::HardClip(i) => (Cigar::HardClip(pos), Cigar::HardClip(*i-pos)),
            Cigar::Pad(i) => (Cigar::Pad(pos), Cigar::Pad(*i-pos)),
            Cigar::Equal(i) => (Cigar::Equal(pos), Cigar::Equal(*i-pos)),
            Cigar::Diff(i) => (Cigar::Diff(pos), Cigar::Diff(*i-pos)),
        }
    }

    /// Returns chopped cigar for chunk size and number of ref bases consumed
    fn chop_cigar(&self, rev_cigar: &mut CigarString) -> (CigarString, u32) {
        let mut num_query_bases_consumed = 0;
        let mut num_ref_bases_consumed = 0;
        let mut cigar_elements = Vec::new();
        while num_query_bases_consumed < self.chunk_size {
            if let Some(c) = rev_cigar.pop() {
                let cigar_pos = (self.chunk_size - num_query_bases_consumed) as u32;
                if Self::consumes_query(&c) {
                    if c.len() <= cigar_pos {
                        cigar_elements.push(c);
                        num_ref_bases_consumed += Self::get_num_ref_bases_consumed(&c);
                    } else {
                        let (c1, c2) = Self::split_cigar(&c, cigar_pos);
                        cigar_elements.push(c1);
                        num_ref_bases_consumed += Self::get_num_ref_bases_consumed(&c1);
                        rev_cigar.push(c2);
                    }
                    num_query_bases_consumed += c.len() as usize;
                } else {
                    cigar_elements.push(c);
                    num_ref_bases_consumed += Self::get_num_ref_bases_consumed(&c);
                }

            } else {
                break;
            }
        }

        (CigarString(cigar_elements), num_ref_bases_consumed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(qname: &str, seq: &str, base_quals: &str, cigar: &CigarString, pos: i64) -> Record {
        let mut rec = Record::default();
        rec.set(qname.as_bytes(), Some(cigar), seq.as_bytes(), base_quals.as_bytes());
        rec.set_pos(pos);
        rec.set_tid(1);
        rec.set_mapq(60);
        rec.set_flags(0);
        rec
    }

    #[test]
    fn simple_test() {
        let chopper_no_edges = AlignmentChopper::new(5, false, false, None);
        let chopper_with_edges = AlignmentChopper::new(5, true, false, None);

        let cigar = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(2), Cigar::Ins(4), Cigar::SoftClip(3)]);
        let rec = make_record("test", "AGTCGATGCATGC", "?!/??50(?/321", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Match(1), Cigar::Ins(4)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 110);

        let cigar3 = CigarString(vec![Cigar::SoftClip(3)]);
        let rec3 = make_record("test-2", "TGC", "321", &cigar3, 111);

        assert_eq!(chopper_no_edges.chop_read(&rec), vec![rec1.clone(), rec2.clone()]);
        assert_eq!(chopper_with_edges.chop_read(&rec), vec![rec1, rec2, rec3]);
    }

    #[test]
    fn test_pos_with_starting_softclip() {
        let chopper_with_edges = AlignmentChopper::new(5, true, false, None);

        let cigar = CigarString(vec![Cigar::SoftClip(4), Cigar::Equal(1), Cigar::Del(4), Cigar::Match(2), Cigar::Ins(4), Cigar::SoftClip(3)]);
        let rec = make_record("test", "AGTCGATGCATGCA", "?!/??50(?/3210", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::SoftClip(4), Cigar::Equal(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Del(4), Cigar::Match(2), Cigar::Ins(3)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 101);

        let cigar3 = CigarString(vec![Cigar::Ins(1), Cigar::SoftClip(3)]);
        let rec3 = make_record("test-2", "TGCA", "3210", &cigar3, 107);
    }

    #[test]
    fn skip_softclips_test() {
        let chopper_skip_softclips_no_edges = AlignmentChopper::new(5, false, true, None);
        let chopper_skip_softclips_with_edges = AlignmentChopper::new(5, true, true, None);

        let cigar = CigarString(vec![Cigar::SoftClip(1), Cigar::Match(4), Cigar::Del(5), Cigar::Match(2), Cigar::Ins(4), Cigar::Equal(1), Cigar::SoftClip(3)]);
        let rec = make_record("test", "CAGTCGATGCATGCG", "??!/??50(?/3210", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Match(1), Cigar::Ins(4)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 110);

        let cigar3 = CigarString(vec![Cigar::Equal(1)]);
        let rec3 = make_record("test-2", "T", "3", &cigar3, 111);

        assert_eq!(chopper_skip_softclips_no_edges.chop_read(&rec), vec![rec1.clone(), rec2.clone()]);
        assert_eq!(chopper_skip_softclips_with_edges.chop_read(&rec), vec![rec1, rec2, rec3]);
    }
}
