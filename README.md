# spadas

## Compile

mvn clean package

## Datasets

put your dataset under the dataset/argo

## run experiments

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.dss.similarity.Framework ./dataset/argo/#your_dataset#

## dataset search
the argo dataset loading: spadas.dss.similarity.Framework.readDatalake()

the search algorithms: spadas.dss.similarity.search

the index structure based on KD-tree: indexStruct.trajectory.clustering.kmeans.indexNode
a

Technical Report
============
https://arxiv.org/abs/2412.04805


Citation
---------
* If you use our code for research work, please cite our paper below:

```
@inproceedings{yang2025,

      title={A Unified Approach for Multi-Granularity Search over Spatial Datasets},

      author={Wenzhe Yang and Sheng Wang and Shixun Huang and Hao Liu and Yuan Sun and Juliana Freire and Zhiyong Peng},

      year={2025},

      eprint={2412.04805},

      archivePrefix={arXiv},

      primaryClass={cs.DB},

      url={https://arxiv.org/abs/2412.04805},

}

