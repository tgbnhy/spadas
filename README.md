# spadas

## Compile

mvn clean package

## Datasets

put your dataset under the dataset/

## run experiments

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar Framework ./dataset/argoverse/data

## dataset search
the argo dataset loading: Framework.readDatalake()

the search algorithms: core/spadas/dss/similarity/Framework

the index structure based on KD-tree: tree/trajectory/clustering/kmeans/indexNode
