# spadas

## Compile

mvn clean package

## Datasets

put your dataset under the dataset/argo

## run experiments

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.dss.similarity.Framework ./dataset/argo/#your_dataset#

## dataset search
the argo dataset loading: edu.nyu.dss.similarity.Framework.readDatalake()

the search algorithms: edu.nyu.dss.similarity.search

the index structure based on KD-tree: au.edu.rmit.trajectory.clustering.kmeans.indexNode
