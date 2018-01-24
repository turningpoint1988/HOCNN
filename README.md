# HOCNN
High-order convolutional neural network for predicting DNA-protein binding sites
This work is based on [caffe-cnn](https://github.com/gifford-lab/caffe-cnn/), and please refer to it for more details.

## Prerequisite
[Caffe and pycaffe](http://caffe.berkeleyvision.org/installation.html)

## Data preparation
	+ Usage example:
	
		```
		python HOCNN.py example/train.tsv example/train_target.tsv example/data/train.h5 -k 2
		```
  + Type the following for details on other optional arguments:
	
		```
		python HOCNN.py -h
		```
