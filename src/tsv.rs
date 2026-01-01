use anyhow::Result;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// A single data point in a TSV track
#[derive(Debug, Clone)]
pub struct TsvPoint
{
    pub chromosome: String,
    pub position: usize,
    pub signal: f64,
}

/// Track for one chromosome
#[derive(Clone)]
pub struct TsvChromosomeTrack
{
    pub points: Vec<TsvPoint>,
    pub min_signal: f64,
    pub max_signal: f64,
}

impl TsvChromosomeTrack
{
    pub fn new() -> Self
    {
        Self {
            points: Vec::new(),
            min_signal: f64::INFINITY,
            max_signal: f64::NEG_INFINITY,
        }
    }

    /// Add a point and update min/max
    pub fn add_point(&mut self, point: TsvPoint)
    {
        if point.signal < self.min_signal
        {
            self.min_signal = point.signal;
        }
        if point.signal > self.max_signal
        {
            self.max_signal = point.signal;
        }
        self.points.push(point);
    }

    /// Sort points by position
    pub fn sort(&mut self)
    {
        self.points.sort_by_key(|p| p.position);
    }

    /// Get points in a region
    pub fn get_points_in_region(&self, start: usize, end: usize) -> Vec<&TsvPoint>
    {
        self.points
            .iter()
            .filter(|p| p.position >= start && p.position <= end)
            .collect()
    }
}

/// Container for all TSV tracks
#[derive(Clone)]
pub struct TsvData
{
    pub tracks: HashMap<String, TsvChromosomeTrack>,
    pub track_name: String,
    pub global_min: f64,
    pub global_max: f64,
}

impl TsvData
{
    /// Create new TSV data container
    pub fn new(track_name: String) -> Self
    {
        Self {
            tracks: HashMap::new(),
            track_name,
            global_min: f64::INFINITY,
            global_max: f64::NEG_INFINITY,
        }
    }

    /// Parse TSV from file (native only)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &str) -> Result<Self>
    {
        use std::fs::File;
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let track_name = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("Custom Track")
            .to_string();

        Self::parse_tsv(reader, track_name)
    }

    /// Parse TSV from bytes
    pub fn from_bytes(data: Vec<u8>, track_name: String) -> Result<Self>
    {
        let reader = BufReader::new(&data[..]);
        Self::parse_tsv(reader, track_name)
    }

    /// Parse TSV data
    fn parse_tsv<R: BufRead>(reader: R, track_name: String) -> Result<Self>
    {
        let mut tsv_data = Self::new(track_name);

        for (line_num, line_result) in reader.lines().enumerate()
        {
            let line = line_result?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#')
            {
                continue;
            }

            // Parse tab-separated fields
            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() < 3
            {
                eprintln!(
                    "Warning: Line {} has fewer than 3 fields, skipping: {}",
                    line_num + 1,
                    line
                );
                continue;
            }

            let chromosome = fields[0].to_string();

            let position = match fields[1].parse::<usize>()
            {
                Ok(pos) => pos,
                Err(e) =>
                {
                    eprintln!(
                        "Warning: Line {} has invalid position '{}': {}",
                        line_num + 1,
                        fields[1],
                        e
                    );
                    continue;
                }
            };

            let signal = match fields[2].parse::<f64>()
            {
                Ok(sig) => sig,
                Err(e) =>
                {
                    eprintln!(
                        "Warning: Line {} has invalid signal '{}': {}",
                        line_num + 1,
                        fields[2],
                        e
                    );
                    continue;
                }
            };

            // Update global min/max
            if signal < tsv_data.global_min
            {
                tsv_data.global_min = signal;
            }
            if signal > tsv_data.global_max
            {
                tsv_data.global_max = signal;
            }

            // Add point to appropriate chromosome track
            let track = tsv_data
                .tracks
                .entry(chromosome.clone())
                .or_insert_with(TsvChromosomeTrack::new);

            track.add_point(TsvPoint {
                chromosome,
                position,
                signal,
            });
        }

        // Sort all tracks by position
        for track in tsv_data.tracks.values_mut()
        {
            track.sort();
        }

        Ok(tsv_data)
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_parse_tsv()
    {
        let data = b"chr1\t100\t1.5\nchr1\t200\t2.3\nchr2\t150\t0.8\n";
        let tsv = TsvData::from_bytes(data.to_vec(), "test".to_string()).unwrap();

        assert_eq!(tsv.tracks.len(), 2);
        assert!(tsv.tracks.contains_key("chr1"));
        assert!(tsv.tracks.contains_key("chr2"));

        let chr1_track = &tsv.tracks["chr1"];
        assert_eq!(chr1_track.points.len(), 2);
        assert_eq!(chr1_track.points[0].position, 100);
        assert_eq!(chr1_track.points[0].signal, 1.5);
        assert_eq!(chr1_track.points[1].position, 200);
        assert_eq!(chr1_track.points[1].signal, 2.3);

        assert_eq!(tsv.global_min, 0.8);
        assert_eq!(tsv.global_max, 2.3);
    }

    #[test]
    fn test_get_points_in_region()
    {
        let mut track = TsvChromosomeTrack::new();
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            position: 100,
            signal: 1.0,
        });
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            position: 200,
            signal: 2.0,
        });
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            position: 300,
            signal: 3.0,
        });

        let points = track.get_points_in_region(150, 250);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].position, 200);
    }
}
