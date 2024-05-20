use std::cmp::min;
use rust_htslib::bam::{Record};
use rust_htslib::bam::record::{CigarString, Cigar, Aux};

#[derive(Debug)]
pub struct AlignmentChopper {
    chunk_size: u32,
    min_length: u32,
    skip_clipped_bases: bool,
    read_group: Option<String>,
    rec_pieces_buffer: Vec<Record>,
    record_slice_meta_buffer: RecordSliceMetaBuffer,
}

#[derive(Debug)]
struct SplitCigarBuf {
    left_c: Cigar,
    right_c: Option<Cigar>,
    query_offset: u32,
    ref_offset: i64,
}

impl SplitCigarBuf {
    fn new(left_c: Cigar, right_c: Option<Cigar>, query_offset: u32, ref_offset: i64) -> Self {
        Self {
            left_c,
            right_c,
            query_offset,
            ref_offset
        }
    }
}

// A struct to hold metadata about current record slicing process
#[derive(Debug)]
struct RecordSliceMetaBuffer {
    global_ref_offset: i64,
    global_query_offset: usize,
    cigar_string: CigarString,
}

impl RecordSliceMetaBuffer {
    fn new() -> Self {
        Self {
            global_ref_offset: 0,
            global_query_offset: 0,
            cigar_string: CigarString(Vec::new()),
        }
    }

    fn reset(&mut self) {
        self.global_ref_offset = 0;
        self.global_query_offset = 0;
        self.cigar_string.clear();
    }
}

impl AlignmentChopper {
    pub fn new(chunk_size: u32, min_length: u32, skip_clipped_bases: bool, read_group: Option<String>) -> Self {
        Self {
            chunk_size,
            min_length,
            skip_clipped_bases,
            read_group,
            rec_pieces_buffer: Vec::new(),
            record_slice_meta_buffer: RecordSliceMetaBuffer::new()
        }
    }

    fn reset(&mut self) {
        // Reset internal buffers for new Record
        self.rec_pieces_buffer.clear();
        self.record_slice_meta_buffer.reset();
    }

    fn add_chunk_record(&mut self, original_rec: &Record, local_query_consumed: usize) {
        let mut new_rec = Record::default();

        // Get seq and qual slices
        let query_offset = self.record_slice_meta_buffer.global_query_offset;
        let chunk_num = self.rec_pieces_buffer.len();
        let slice_end = min(original_rec.seq_len(), query_offset + local_query_consumed);

        let new_seq = &original_rec.seq().as_bytes()[query_offset..slice_end];
        let new_qual = &original_rec.qual()[query_offset..slice_end];

        // Update name for chunk
        let new_qname = &[original_rec.qname(), b"-", chunk_num.to_string().as_bytes()].concat();

        // These are changed based on the particular slice
        new_rec.set(new_qname, Some(&self.record_slice_meta_buffer.cigar_string), new_seq, new_qual);
        new_rec.set_pos(original_rec.pos() + self.record_slice_meta_buffer.global_ref_offset);

        // Following are unchanged
        new_rec.set_flags(original_rec.flags());
        new_rec.set_tid(original_rec.tid());
        new_rec.set_mapq(original_rec.mapq());
        new_rec.set_mtid(original_rec.mtid());
        new_rec.set_mpos(original_rec.mpos());
        new_rec.set_insert_size(original_rec.insert_size());

        // All aux data other than RG is lost
        if let Some(rg) = &self.read_group {
            if let Ok(_a) = new_rec.aux(b"RG") {
                new_rec.remove_aux(b"RG").expect(&format!("Could not remove RG from: {} - {}", &new_rec.tid(), &new_rec.pos()));
            }
            new_rec.push_aux(b"RG", Aux::String(rg)).expect(&format!("Unable to push RG string at: {} - {}", &new_rec.tid(), &new_rec.pos()));
        }

        self.rec_pieces_buffer.push(new_rec);
    }

    fn consume_cigar(c: &Cigar, amount: u32) -> SplitCigarBuf {
        match c {
            Cigar::Match(x) => {
                let (left_c, right_c, consumed) = if amount < *x {
                    (Cigar::Match(amount), Some(Cigar::Match(*x - amount)), amount)
                } else {
                    (Cigar::Match(*x), None, *x)
                };
                SplitCigarBuf::new(left_c, right_c, consumed, consumed as i64)
            },
            Cigar::Ins(x) => {
                let (left_c, right_c, consumed) = if amount < *x {
                    (Cigar::Ins(amount), Some(Cigar::Ins(*x - amount)), amount)
                } else {
                    (Cigar::Ins(*x), None, *x)
                };
                SplitCigarBuf::new(left_c, right_c, consumed, 0i64)
            },
            Cigar::Del(x) => {
                SplitCigarBuf::new(Cigar::Del(*x), None, 0, *x as i64)
            },
            Cigar::RefSkip(x) => {
                SplitCigarBuf::new(Cigar::RefSkip(*x), None, 0, *x as i64)
            },
            Cigar::SoftClip(x) => {
                let (left_c, right_c, consumed) = if amount < *x {
                    (Cigar::SoftClip(amount), Some(Cigar::SoftClip(*x - amount)), amount)
                } else {
                    (Cigar::SoftClip(*x), None, *x)
                };
                SplitCigarBuf::new(left_c, right_c, consumed, 0i64)
            },
            Cigar::HardClip(x) => {
                SplitCigarBuf::new(Cigar::HardClip(*x), None, 0, 0i64)
            },
            Cigar::Pad(x) => {
                SplitCigarBuf::new(Cigar::Pad(*x), None, 0, 0i64)
            },
            Cigar::Equal(x) => {
                let (left_c, right_c, consumed) = if amount < *x {
                    (Cigar::Equal(amount), Some(Cigar::Equal(*x - amount)), amount)
                } else {
                    (Cigar::Equal(*x), None, *x)
                };
                SplitCigarBuf::new(left_c, right_c, consumed, consumed as i64)
            }
            Cigar::Diff(x) => {
                let (left_c, right_c, consumed) = if amount < *x {
                    (Cigar::Diff(amount), Some(Cigar::Diff(*x - amount)), amount)
                } else {
                    (Cigar::Diff(*x), None, *x)
                };
                SplitCigarBuf::new(left_c, right_c, consumed, consumed as i64)
            }
        }
    }

    pub fn chop_read(&mut self, rec: &Record) -> &Vec<Record> {
        self.reset();  // Clear internal buffers

        let mut local_ref_consumed = 0;
        let mut local_query_consumed = 0;

        let mut cigar_consumption;

        let mut current_cigar = rec.cigar().take();

        // Handle trailing clipped bases
        if self.skip_clipped_bases {
            if rec.cigar().trailing_hardclips() > 0 {
                current_cigar.pop();
            }
            let trailing_softclips = rec.cigar().trailing_softclips() as usize;
            if trailing_softclips > 0 {
                current_cigar.pop();
            }
        }

        let mut cigar_iter = current_cigar.into_iter();

        // Handle starting clipped bases
        if self.skip_clipped_bases {
            if rec.cigar().leading_hardclips() > 0 {
                cigar_iter.next();
            }
            let leading_softclips = rec.cigar().leading_softclips() as usize;
            if leading_softclips > 0 {
                cigar_iter.next();
                self.record_slice_meta_buffer.global_query_offset += leading_softclips;
            }
        }

        for c in cigar_iter {
            cigar_consumption = Self::consume_cigar(c, self.chunk_size - local_query_consumed);
            self.record_slice_meta_buffer.cigar_string.push(cigar_consumption.left_c);
            local_ref_consumed += cigar_consumption.ref_offset;
            local_query_consumed += cigar_consumption.query_offset;

            if cigar_consumption.right_c.is_none() {
                // Fully consumed cigar token
                if local_query_consumed == self.chunk_size {
                    // Add record if filled chunk_size
                    self.add_chunk_record(rec, local_query_consumed as usize);

                    // Update global offsets after adding records
                    self.record_slice_meta_buffer.global_ref_offset += local_ref_consumed;
                    self.record_slice_meta_buffer.global_query_offset += local_query_consumed as usize;

                    // Restart new consumption cycle
                    self.record_slice_meta_buffer.cigar_string.clear();
                    local_ref_consumed = 0;
                    local_query_consumed = 0;
                }
            } else {
                // Finish consuming any Cigar in the buffer from previous iteration
                while let Some(c_buf) = cigar_consumption.right_c {
                    // Partially consumed cigar, so must be time to write new record chunk
                    self.add_chunk_record(rec, local_query_consumed as usize);

                    // Update global offsets after adding records
                    self.record_slice_meta_buffer.global_ref_offset += local_ref_consumed;
                    self.record_slice_meta_buffer.global_query_offset += local_query_consumed as usize;

                    // Restart new consumption cycle
                    self.record_slice_meta_buffer.cigar_string.clear();
                    local_ref_consumed = 0;
                    local_query_consumed = 0;

                    cigar_consumption = Self::consume_cigar(&c_buf, self.chunk_size - local_query_consumed);
                    self.record_slice_meta_buffer.cigar_string.push(cigar_consumption.left_c);
                    local_ref_consumed += cigar_consumption.ref_offset;
                    local_query_consumed += cigar_consumption.query_offset;
                }
            }
        }

        // Handle min length requirement for last chunk
        if local_query_consumed >= self.min_length {
            self.add_chunk_record(rec, local_query_consumed as usize);
        }

        &self.rec_pieces_buffer
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
        let mut chopper_no_edges = AlignmentChopper::new(5, 5, false, None);
        let mut chopper_with_edges = AlignmentChopper::new(5, 0, false, None);

        let cigar = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(2), Cigar::Ins(4), Cigar::SoftClip(3)]);
        let rec = make_record("test", "AGTCGATGCATGC", "?!/??50(?/321", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Match(1), Cigar::Ins(4)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 110);

        let cigar3 = CigarString(vec![Cigar::SoftClip(3)]);
        let rec3 = make_record("test-2", "TGC", "321", &cigar3, 111);

        assert_eq!(chopper_no_edges.chop_read(&rec), &vec![rec1.clone(), rec2.clone()]);
        assert_eq!(chopper_with_edges.chop_read(&rec), &vec![rec1, rec2, rec3]);
    }

    #[test]
    fn test_pos_with_starting_softclip() {
        let mut chopper_with_edges = AlignmentChopper::new(5, 0, false, None);

        let cigar = CigarString(vec![Cigar::SoftClip(4), Cigar::Equal(1), Cigar::Del(4), Cigar::Match(2), Cigar::Ins(4), Cigar::SoftClip(3)]);
        let rec = make_record("test", "AGTCGATGCATGCA", "?!/??50(?/3210", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::SoftClip(4), Cigar::Equal(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Del(4), Cigar::Match(2), Cigar::Ins(3)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 101);

        let cigar3 = CigarString(vec![Cigar::Ins(1), Cigar::SoftClip(3)]);
        let rec3 = make_record("test-2", "TGCA", "3210", &cigar3, 107);

        assert_eq!(chopper_with_edges.chop_read(&rec), &vec![rec1, rec2, rec3]);
    }

    #[test]
    fn skip_softclips_test() {
        let mut chopper_skip_softclips_no_edges = AlignmentChopper::new(5, 5, true, None);
        let mut chopper_skip_softclips_with_edges = AlignmentChopper::new(5, 0, true, None);

        let cigar = CigarString(vec![Cigar::SoftClip(1), Cigar::Match(4), Cigar::Del(5), Cigar::Match(2), Cigar::Ins(4), Cigar::Equal(1), Cigar::SoftClip(3)]);
        let rec = make_record("test", "CAGTCGATGCATGCG", "??!/??50(?/3210", &cigar, 100);

        let cigar1 = CigarString(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(1)]);
        let rec1 = make_record("test-0", "AGTCG", "?!/??", &cigar1, 100);

        let cigar2 = CigarString(vec![Cigar::Match(1), Cigar::Ins(4)]);
        let rec2 = make_record("test-1", "ATGCA", "50(?/", &cigar2, 110);

        let cigar3 = CigarString(vec![Cigar::Equal(1)]);
        let rec3 = make_record("test-2", "T", "3", &cigar3, 111);

        assert_eq!(chopper_skip_softclips_no_edges.chop_read(&rec), &vec![rec1.clone(), rec2.clone()]);
        assert_eq!(chopper_skip_softclips_with_edges.chop_read(&rec), &vec![rec1, rec2, rec3]);
    }

    #[test]
    fn large_clips_test() {

    }
}
