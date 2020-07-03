#!/usr/local/bin/python3
import sys,json,datetime,os
# -crf sets the compression 0=lossless 51=max 23=default (see https://trac.ffmpeg.org/wiki/Encode/H.264)
#
# hqdn3d is the denoise filter but there are many others, see 
#
# -preset ultrafast|superfast|veryfast|faster|fast|medium – default|slow|slower|veryslow|placebo
# veryfast is about twice faster then mediim (default)
# the slower is the preset, the better is the quality of the video
#
# si potrebbe usare anche la seguente opzione
# -tune     film – use for high quality movie content; lowers deblocking
#           animation – good for cartoons; uses higher deblocking and more reference frames
#           grain – preserves the grain structure in old, grainy film material
#           stillimage – good for slideshow-like content
#           fastdecode – allows faster decoding by disabling certain filters
#           zerolatency – good for fast encoding and low-latency streaming
#           psnr – ignore this as it is only used for codec development
#           ssim – ignore this as it is only used for codec development 

prog='ffmpeg'

#helper functions
def print_error():
    print('syntax: extract_data.py [--compression|-c <0-51 0=lossless>')
    print(' -s|--speed <ultrafast|superfast|veryfast|faster|fast|medium|slow|slower|veryslow|placebo >')
    print(' -tune  <film|animation|grain|stillimage|fastdecode|zerolatency|psnr|ssim>')
    print(' -co|--codec <x264|x265> --filter|-f <hqdn|ata|nlm|vag|bm3d> --deflicker|-d] <input_file> <output_file default=out.mkv>')

args = sys.argv
#print(args)
if len(args) < 2:
    print_error()
    quit()
itargs=iter(args)

infile= ''
outfile= 'out.mkv'
nfilter = 'hqdn3d'
speed = 'medium'
codec= 'x264'
defl=False
compress='26' #default is 23 ma this way file has a size equal to original nikon MOV
del(args[0])
for a in itargs:
    if a == '--help' or a  == '-h':
        print_error()
        quit()
    elif a == '-c' or a == '--compression':
        compress=next(itargs)
    elif a == '-s' or a == '--speed':
        speed=next(itargs)
    elif a == '-co' or a == '--codec':
        codec=next(itargs)
    elif a == '-f' or a == '--filter':
        nfilter = next(itargs)
    elif a == '-d' or a == '--deflicker':
        defl = True
    else:
        print('a=', a)
        infile = a
        outfile= next(itargs)
#hqdn3dpars='=4.0:3.0:6.0:4.5' # this values are the default ones
hqdn3dpars='' # default values

#use hqdn (this should be the best one!) or ata which are faster
#vaguedenoiser is good bu slower
#nlmeans and bm3d are way too slow...
# for further details see https://ffmpeg.org/ffmpeg-filters.html#commands
if nfilter == 'hqdn':
    denstr='hqdn3d'+hqdn3dpars
elif nfilter == 'ata':
    denstr='atadenoise'
elif nfilter == 'nlm':
    denstr = 'nlmeans' # slow but very good
elif nfilter == 'vag':
    denstr = 'vaguedenoiser'
elif nfilter == 'bm':
    denstr = 'bm3d'
else:
    denstr='hqdn3d'+hqdn3dpars

if codec=='x264':
    codstr='libx264'
elif compress=='x265':
    codstr='libx265'
else:
    codstr='libx264'
exestr=prog + ' -i ' + infile + ' -c:v ' + codstr + ' -c:a copy' + ' -crf ' + compress + ' -preset ' + speed + ' -filter:v ' + denstr
if defl==True:
    exestr = exestr + ' -filter:v deflicker'
exestr = exestr + ' ' + outfile
print(exestr)
os.system(exestr)
#ffmpeg -i infile -c:v libx264 -c:a copy -crf 0 -preset ultrafast -filter:v hqdn3d=4.0:3.0:6.0:4.5 $2
