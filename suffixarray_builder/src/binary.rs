use std::cmp::min;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};

const ONE_GIB: usize = 2usize.pow(30);

pub trait Serializable {
    fn serialize(&self) -> Vec<u8>;
}

impl Serializable for [i64] {
    fn serialize(&self) -> Vec<u8> {
        let mut res = vec![];
        self.iter().for_each(|entry|
            res.extend_from_slice(&entry.to_le_bytes())
        );
        res
    }
}

fn deserialize_sa(data: &[u8]) -> Vec<i64> {
    let mut res = vec![];
    if data.len() % 8 != 0 {
        panic!("Serialized data is not a multiple of 8 bytes long!")
    }
    for start in (0..data.len()).step_by(8) {
        res.push(i64::from_le_bytes(data[start..start + 8].try_into().unwrap()));
    }
    res
}

// from: https://gist.github.com/taylorsmithgg/ba7b070c0964aa8b86d311ab6f8f5508
// https://dev.to/oliverjumpertz/how-to-write-files-in-rust-m06?comments_sort=top
pub fn write_binary(sample_rate: u8, suffix_array: &[i64], name: &str) -> Result<(), std::io::Error> {
    // TODO: how to store the uniprot protein data? store in separate files, or just assume we re-read the complete tsv
    //  we could also use the first x bytes to store the uniprot version that this SA was built for
    // create the file
    let mut f = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true) // if the file already exists, empty the file
        .open(name.to_owned() + "_sa.bin")?;
    f.write_all(&[sample_rate])?; // write the sample rate as the first byte

    // write 1 GiB at a time, to minimize extra used memory since we need to translate i64 to [u8; 8]
    let sa_len = suffix_array.len();
    for start_index in (0..sa_len).step_by(ONE_GIB/8) {
        let end_index = min(start_index + ONE_GIB/8, sa_len);
        f.write_all(&suffix_array[start_index..end_index].serialize())?;
    }

    Ok(())
}

pub fn load_binary(name: &str) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    // read the SA and deserialize it into a vec of i64
    let sa_file = File::open(name)?;
    let (sample_rate, sa) = read_sa_file(&sa_file)?;
    Ok((sample_rate, sa))
}

fn read_sa_file(mut file: &File) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    let mut sample_rate_buffer = [0_u8; 1];
    file.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    // this buffer is 1GiB big
    let mut sa = vec![];
    loop {
        let mut buffer = vec![];
        // use take in combination with read_to_end to ensure that the buffer will be completely filled (except when the file is smaller than the buffer)
        let count = file.take(ONE_GIB as u64).read_to_end(&mut buffer)?;
        if count == 0 {
            break;
        }
        sa.extend_from_slice(&deserialize_sa(&buffer[..count]));
    }

    Ok((sample_rate, sa))
}


#[cfg(test)]
mod tests {
    use crate::binary::{deserialize_sa, Serializable};

    #[test]
    fn test_serialize_deserialize() {
        let data: Vec<i64> = vec![5, 2165487362, -12315135];
        let serialized = data.serialize();
        let deserialized = deserialize_sa(serialized.as_ref());
        assert_eq!(data, deserialized);
    }

    #[test]
    fn test_serialize_deserialize_empty() {
        let data: Vec<i64> = vec![];
        let serialized = data.serialize();
        let deserialized = deserialize_sa(serialized.as_ref());
        assert_eq!(data, deserialized);
    }
}