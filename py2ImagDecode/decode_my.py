#coding=utf-8
import json,re,sys,os,logging,operator,numpy,random
from reedsolo import RSCodec
from collections import defaultdict
from shutil import copyfile
from utils.robust_solition import PRNG
from utils.droplet import Droplet
import utils.file_process as fp
from tqdm import tqdm
import utils.Colorer
import json,re,sys,os,logging,operator,numpy,random
from time import sleep
import hashlib

logging.basicConfig(level=logging.DEBUG)
sys.setrecursionlimit(10000000)
rs = 5
rscode = RSCodec(rs)
max_hamming = 100
header_size = 4
chunk_num = 1494
chunks = [None] * chunk_num
prng = PRNG(K = chunk_num, delta = 0.05, c = 0.1, np = False)
droplets = set()
done_segments = set()
chunk_to_droplets = defaultdict(set)

def get_payload(dna):
    # print(dna)
    # convert a string like ACTCA to an array of ints like [10, 2, 4]
    num = dna.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3') # 一个碱基编码2个bit，一个信息byte是8个bit
    s = ''.join('{0:02b}'.format(int(num[t])) for t in range(0, len(num), 1))  ##转成二进制
    # print(s)
    data = [int(s[t:t + 8], 2) for t in range(0, len(s), 8)] #每四个为一组的二进制变成一个10进制数，相当于每8个bit为一个信息点，将其转化为一个10进制数
    # print(data)
    ## 有错误,如果有就纠正
    try:
        # print(rscode.decode(data)[0])
        data_corrected = list(rscode.decode(data)[0])
        # print(data_corrected)
    except:
        return -1, None, None  # could not correct the code
    data_again = list(rscode.encode(data_corrected))
    ## 如果经过RS修正后解码结果与原数据差距过大
    if numpy.count_nonzero(list(data[i] != data_again[i] for i in range(0,len(data_corrected),1))) > max_hamming:
        ## 太多错误了,直接返回跳出
        return -1, None, None
    # get the payload
    seed_array = data_corrected[:header_size]
    seed = sum([int(x) * 256 ** i for i, x in enumerate(seed_array[::-1])])
    # print(seed)
    payload = data_corrected[header_size:]
    prng.set_seed(seed)
    ix_samples = prng.get_src_blocks_wrap()[1]
    # print(ix_samples)
    return payload, seed, ix_samples

def addDroplet(droplets, chunk_to_droplets, droplet, done_segments):
    droplets.add(droplet)
    for chunk_num in droplet.num_chunks:
        chunk_to_droplets[chunk_num].add(droplet)
    updateEntry(droplets,chunk_to_droplets,droplet,done_segments)

def updateEntry(droplets, chunk_to_droplets, droplet, done_segments):

    for chunk_num in (droplet.num_chunks & done_segments):
        droplet.data = list(map(operator.xor, droplet.data, chunks[chunk_num]))
        droplet.num_chunks.remove(chunk_num)
        chunk_to_droplets[chunk_num].discard(droplet)

    if len(droplet.num_chunks) == 1:
        lone_chunk = droplet.num_chunks.pop()

        chunks[lone_chunk] = droplet.data

        done_segments.add(lone_chunk)
        droplets.discard(droplet)
        chunk_to_droplets[lone_chunk].discard(droplet)

        for other_droplet in chunk_to_droplets[lone_chunk].copy():
            updateEntry(droplets, chunk_to_droplets, other_droplet, done_segments)

def main():
    file = '50-SF.txt'
    outfile = 'output.jpg'
    original_file =  '50-SF.jpg'
    count_sleep = 1
    with open(file, 'r') as f:
        dna = f.readline().rstrip('\n')
        while (dna):
            payload, seed, ix_samples = get_payload(dna)
            droplet = Droplet(payload, seed, ix_samples)
            addDroplet(droplets, chunk_to_droplets, droplet, done_segments)
            if (chunk_num - len(done_segments) <= 0):
                break
            dna = f.readline().rstrip()
        outstring = ''
        logging.info("We are restoring the picture now!")
        for x in tqdm(chunks):
            outstring += ''.join(map(chr, x))
            count_sleep +=1
            if count_sleep % 9 == 0:
                sleep(0.001)
        with open(outfile, 'rb+') as f:
            f.write(outstring)
            md_decode = hashlib.md5()
            md_decode.update(f.read())
            out_decode = md_decode.hexdigest()
            logging.info("The MD5 for decoded picture is %s", fp.get_md5(outstring))
        with open(original_file, 'rb') as f:
            all_lines = f.read()
            md = hashlib.md5()
            md.update(f.read())
            out = md.hexdigest()
            logging.info("The MD5 for decoded picture in code offered by teacher is %s",fp.get_md5(all_lines))

if __name__=='__main__':
    main()

