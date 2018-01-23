import argparse,pwd,os,numpy as np,h5py,sys,math
from os.path import splitext,exists,dirname
from os import makedirs
from itertools import izip

def outputHDF5(data,label,filename,labelname,dataname):
    print 'data shape: ',data.shape
    comp_kwargs = {'compression': 'gzip', 'compression_opts': 1}
    label = [[x.astype(np.float32)] for x in label]
    with h5py.File(filename, 'w') as f:
    	f.create_dataset(dataname, data=data, **comp_kwargs)
    	f.create_dataset(labelname, data=label, **comp_kwargs)

def seq2feature(data,mapper,label,out_filename,labelname,dataname):
	out = []
	for seq in data:
		mat = embed(seq,mapper)
		result = mat.transpose()
		result1 = [ [a] for a in result]
		out.append(result1)

	outputHDF5(np.asarray(out),label,out_filename,labelname,dataname)

def embed(seq,mapper):
	mat = []
	for element in seq:
	   if mapper.has_key(element):
	      mat.append(mapper.get(element))
	   else:
	      print >> sys.stderr,"invalid character!";sys.exit(1)
	return np.asarray(mat)

def convert(infile,labelfile,outfile,mapper,batchsize,labelname,dataname,degree):
	with open(infile) as seqfile, open(labelfile) as labelfile:
		cnt = 0
		seqdata = []
		label = []
		batchnum = 0
		for x,y in izip(seqfile,labelfile):
			temp = []
			sequence = x.strip().split()[1]
			length = len(sequence)
			for i in range(length-degree+1):
				multinucleotide = sequence[i:i+degree]
				temp.append(multinucleotide)
			seqdata.append(temp)
			label.append(float(y.strip()))
			cnt = (cnt+1) % batchsize
			if cnt == 0:
				batchnum = batchnum + 1
				seqdata = np.asarray(seqdata)
				label = np.asarray(label)
				t_outfile = outfile + '.batch' + str(batchnum)
				seq2feature(seqdata,mapper,label,t_outfile,labelname,dataname)
				seqdata = []
				label = []
		
		if cnt >0:
			batchnum = batchnum + 1
			seqdata = np.asarray(seqdata)
			label = np.asarray(label)
			t_outfile = outfile + '.batch' + str(batchnum)
			seq2feature(seqdata,mapper,label,t_outfile,labelname,dataname)

	return batchnum


# build a code table
def buildmapper(degree):
	length = degree
	alphabet = ['A','C','G','T']
	mapper = ['']
	while length > 0:
		mapper_len = len(mapper)
		temp = mapper
		for base in range(len(temp)):
			for letter in alphabet:
				mapper.append(temp[base] + letter)
        #delete the original conents
		while mapper_len > 0:
			mapper.pop(0)
			mapper_len -= 1
		
		length -= 1

	code = np.eye(len(mapper), dtype = int)
	encoder = {}
	for i in range(len(mapper)):
		encoder[mapper[i]] = list(code[i,:])
    
	number = int(math.pow(4,degree))
	encoder['N'] = [0]*len(mapper)*number
	return encoder

 
def parse_args():
	parser = argparse.ArgumentParser(description="Convert sequence and target for Caffe")
	user = pwd.getpwuid(os.getuid())[0]

    # Positional (unnamed) arguments:
	parser.add_argument("infile",  type=str, help="Sequence in FASTA/TSV format (with .fa/.fasta or .tsv extension)")
	parser.add_argument("labelfile",  type=str,help="Label of the sequence. One number per line")
	parser.add_argument("outfile",  type=str, help="Output file (example: $MODEL_TOPDIR$/data/train.h5). ")

    # Optional arguments:
	parser.add_argument("-m", "--mapperfile", dest="mapperfile", default="", help="A TSV file mapping each nucleotide to a vector. The first column should be the nucleotide, and the rest denote the vectors.")
	parser.add_argument("-k", "--degree", dest="degree", type=int, default=1, help="The degree of high-order.")
	parser.add_argument("-b", "--batch", dest="batch", type=int, default=5000, help="Batch size for data storage (Defalt:5000)")
	parser.add_argument("-l", "--labelname", dest="labelname",default='label', help="The group name for labels in the HDF5 file")
	parser.add_argument("-d", "--dataname", dest="dataname",default='data', help="The group name for data in the HDF5 file")

	return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    outdir = dirname(args.outfile)
    if not exists(outdir):
        makedirs(outdir)

    if args.mapperfile == "":
        mapper = buildmapper(args.degree)
    else:
        args.mapper = {}
        with open(args.mapperfile,'r') as f:
            for x in f:
                line = x.strip().split()
                word = line[0]
                vec = [float(item) for item in line[1:]]
                args.mapper[word] = vec
                
    batchnum = convert(args.infile,args.labelfile,args.outfile,mapper,args.batch,args.labelname,args.dataname,args.degree)
