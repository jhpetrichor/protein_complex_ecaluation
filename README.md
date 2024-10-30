# protein_complex_ecaluation

## 实现的指标
```rust
pub struct IndexResult {
    precision: f64,     // precision 准确率
    recall: f64,        // recall 召回率
    f_measure: f64,     // f-measure 
    sn: f64,            // sensitivity 灵敏度
    ppv: f64,           // Positive Predictive Value 阳性预测值
    acc: f64            // Accuracy 准确性
}
```

## 配置项
- 用户需要配置PPI，标准蛋白质复合物，预测蛋白质复合物的文件路径
- 将过滤掉小于**min_size**的复合物
- 复合物匹配阈值**threshold**,两个复合物的重叠系数（overlapping score >= threshold, 则认为这两个复合物是匹配的
```rust
pub struct Options {
    pub ppi_path: PathBuf,
    pub ref_complex_path: PathBuf,
    pub pre_complex_path: PathBuf,
    pub min_size: u32,
    pub threshold: f64,
}
```
## example
```rust

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

```
### run
- install Rust: https://www.rust-lang.org/
```shell
cargo run
```
