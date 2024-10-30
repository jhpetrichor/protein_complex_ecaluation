use core::panic;
use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;

use crate::{complex::Complex, options::Options};

/// 存放预测和标准蛋白质复合物，便于分析
pub struct Analysis {
    ref_complex: Vec<Complex>,      // 标准蛋白质复合物
    pre_complex: Vec<Complex>,      // 预测的蛋白质复合物
    os_matrix: Vec<Vec<f64>>,       // 存储复合物之间的匹配指数
    overlap_matrix: Vec<Vec<u32>>,  // 存储复合物之间重叠的蛋白质个数
}

/// 存放结果指标
#[derive(Debug)]
pub struct IndexResult {
    precision: f64,     // precision 准确率
    recall: f64,        // recall 召回率
    f_measure: f64,     // f-measure 
    sn: f64,            // sensitivity 灵敏度
    ppv: f64,           // Positive Predictive Value 阳性预测值
    acc: f64            // Accuracy 准确性
}

impl Analysis {
    pub fn new(options: &Options) -> Self {
        let proteins = read_ppi_proteins_from_file(options.ppi_path.clone());
        let ref_complex = read_complex_from_file(
            options.ref_complex_path.clone(),
            &proteins,
            options.min_size,
        );
        let pre_complex = read_complex_from_file(
            options.pre_complex_path.clone(),
            &proteins,
            options.min_size,
        );
        let (os_matrix, overlap_matrix) = calculate_all_os_overlap(&pre_complex, &ref_complex);

        Self {
            ref_complex,
            pre_complex,
            os_matrix,
            overlap_matrix,
        }
    }

    pub fn calculate_index(&self, threshold: f64) -> IndexResult {
        // 计算准确率
        let mut predicted = 0;
        for row in self.os_matrix.iter() {
            for os in row {
                if os > &threshold {
                    predicted += 1;
                    break;
                }
            }
        }
        let precision = predicted as f64 / self.pre_complex.len() as f64;
        // 计算召回率
        let mut predicted = 0;
        for i in 0..self.ref_complex.len() {
            for j in 0..self.pre_complex.len() {
                if self.os_matrix[j][i] > threshold {
                    predicted += 1;
                    break;
                }
            }
        }
        let recall = predicted as f64 / self.ref_complex.len() as f64;
        // 计算F-measure
        let f_measure = (2.0 * precision * recall) / (precision + recall);

        // 计算Sn
        let sum_count = self.ref_complex.iter().fold(0, |acc, x| acc + x.size());
        let mut sum_max = 0;
        for i in 0..self.ref_complex.len() {
            let mut max = 0;
            for j in 0..self.pre_complex.len() {
                if self.overlap_matrix[j][i] > max {
                    max = self.overlap_matrix[j][i];
                }
            }
            sum_max += max;
        }
        let sn = sum_max as f64 / sum_count as f64;

        // 计算PPV acc
        let sum_count = self.overlap_matrix.iter().fold(0, |acc, row| {
            acc + row.iter().fold(0, |acc1, complex| acc1 + complex)
        });
        let mut sum_max = 0;
        for i in 0..self.pre_complex.len() {
            let mut max = 0;
            for j in 0..self.ref_complex.len() {
                if self.overlap_matrix[i][j] > max {
                    max = self.overlap_matrix[i][j];
                }
            }
            sum_max += max;
        }
        let ppv = sum_max as f64 / sum_count as f64;
        let acc = (sn * ppv).sqrt();

        IndexResult {
            precision,
            recall,
            f_measure,
            sn,
            ppv,
            acc
        }
    }
}

// 从PPI文件中读取蛋白质
fn read_ppi_proteins_from_file(file_path: PathBuf) -> BTreeSet<String> {
    if !file_path.is_file() {
        panic!("ppi_file is empty!")
    }

    let mut proteins = BTreeSet::new();
    let contents = fs::read_to_string(file_path).expect("failed to open ppi file!");
    for line in contents.lines() {
        let pair = line.split_whitespace().collect::<Vec<&str>>();
        proteins.insert(pair[0].to_string());
        proteins.insert(pair[1].to_string());
    }

    proteins
}

// 从文件中读取蛋白质复合物
fn read_complex_from_file(
    file_path: PathBuf,
    proteins: &BTreeSet<String>,
    min_size: u32,
) -> Vec<Complex> {
    let mut complexes = Vec::new();
    let contents = fs::read_to_string(file_path).expect("failed to open complex file!");
    for line in contents.lines() {
        println!("{}", line);
        let complex = line.split_whitespace().collect::<Vec<&str>>();
        if complex.len() < min_size as usize {
            continue;
        }
        let tmp: BTreeSet<String> = complex.iter().map(|p| p.to_string()).collect();
        complexes.push(Complex { proteins: tmp });
    }

    complexes
}

// 计算预测蛋白质复合物和标准蛋白质复合物之间的重叠系数
fn calculate_all_os_overlap(
    pre_complex: &Vec<Complex>,
    ref_complex: &Vec<Complex>,
) -> (Vec<Vec<f64>>, Vec<Vec<u32>>) {
    let mut os_matrix = Vec::with_capacity(pre_complex.len());
    let mut overlap_matrix = Vec::with_capacity(pre_complex.len());
    for _ in 0..ref_complex.len() {
        let row = vec![0.0; ref_complex.len()];
        os_matrix.push(row);
        let row = vec![0; ref_complex.len()];
        overlap_matrix.push(row);
    }

    for (i, p) in pre_complex.iter().enumerate() {
        for (j, r) in ref_complex.iter().enumerate() {
            let (score, size) = p.os(r);
            os_matrix[i][j] = score;
            os_matrix[j][i] = score;
            overlap_matrix[i][j] = size;
            overlap_matrix[j][i] = size;
        }
    }

    (os_matrix, overlap_matrix)
}

#[cfg(test)]
mod tests {
    use crate::options::Options;
    use std::path::PathBuf;

    use super::Analysis;

    #[test]
    fn test_read_complex() {
        let ppi_file = PathBuf::from("./data/ppi.txt");
        let complex_file = PathBuf::from("./data/complex.txt");
        let option = Options {
            ppi_path: ppi_file,
            ref_complex_path: complex_file.clone(),
            pre_complex_path: complex_file,
            min_size: 3,
            threshold: 0.25,
        };
        let analysis = Analysis::new(&option);
        println!("{:?}", analysis.ref_complex);
        println!("{:?}", analysis.pre_complex);
        println!("{:?}", analysis.os_matrix);
        println!("{:?}", analysis.overlap_matrix);

        let index_result = analysis.calculate_index(option.threshold);
        println!("{:?}", index_result);
    }
}
