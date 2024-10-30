use std::path::PathBuf;
use analysis::Analysis;
use options::Options;

mod analysis;
mod complex;
mod options;


fn main() {
    // 初始化配置项
    let option = Options {
        ppi_path: PathBuf::from("path/to/ppi"),
        ref_complex_path: PathBuf::from("path/to/reference"),
        pre_complex_path: PathBuf::from("path/to/predict"),
        min_size: 3,
        threshold: 0.25,
    };

    // 预测结果分析
    let analysis = Analysis::new(&option);
    let index_result = analysis.calculate_index(option.threshold);

    // 输出结果
    println!("{:?}", index_result);
}
