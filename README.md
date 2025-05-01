# SPADAS

## Compile

mvn clean package

## Datasets

put your dataset under the dataset/

| __Dataset__ |__link__|
|-------------|--------|
| MultiOpen   |        |
| T-drive     |https://www.microsoft.com/en-us/research/publication/t-drive-trajectory-data-sample/     |
| Argoverse   |    https://www.argoverse.org/av1.html     |
| ShapeNet    |    https://shapenet.org/    |
| Chicago     |    https://data.cityofchicago.org/Transportation/Taxi-Trips-2013-2023-/wrvz-psew    |
| Proto       |  https://archive.ics.uci.edu/dataset/339/taxi+service+trajectory+prediction+challenge+ecml+pkdd+2015      |

## run experiments

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar Framework ./dataset/argoverse/data

## dataset search
the argo dataset loading: Framework.readDatalake()

the search algorithms: core/spadas/dss/similarity/Framework

the index structure based on KD-tree: tree/trajectory/clustering/kmeans/indexNode

Technical Report
============
https://arxiv.org/abs/2412.04805


Citation
---------
* If you use our code for research work, please cite our paper below:

```
@inproceedings{yang2025,

      title={A Unified Approach for Multi-Granularity Search over Spatial Datasets},

      author={Wenzhe Yang and Sheng Wang and Shixun Huang and Yuyang Liao and Yuan Sun and Juliana Freire and Zhiyong Peng},

      year={2025},

      eprint={2412.04805},

      archivePrefix={arXiv},

      primaryClass={cs.DB},

      url={https://arxiv.org/abs/2412.04805},

}

