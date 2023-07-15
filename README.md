# spadas

## Compile

mvn clean package

## Datasets

put your dataset under the dataset/

## run experiments

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar Framework ./dataset/#your_dataset#

## dataset search
the argo dataset loading: Framework.readDatalake()

the search algorithms: core.spadas.dss.similarity.search

the index structure based on KD-tree: au/core/trajectory/clustering/kmeans/indexNode
