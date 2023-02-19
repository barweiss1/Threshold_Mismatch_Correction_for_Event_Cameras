# Threshold_Mismatch_Correction_for_Event_Cameras
 
This repository contains code regarding the threshold mismatch effect in event cameras. \
The code containts 3 main things:
1. A tool that converts syntheticly generated high temporal resolution videos to event streams.
2. Analysis of a threshold estimation algorithm based on event streams.
3. An implementation and analysis of a novel threshold mismatch correction algorithm. 

The repository also contains a report explaining the work done here inside the presentations folder.

## V2E tool 
**Code/our_v2e/our_v2e_main.m** - contains conversion of video to event streams. \
The tool is based on the tool developed by Yuhuang Hu et al in the paper "v2e: From Video Frames to Realistic DVS Events". \
It was simplified to only include the threshold mismatch effect and eliminate other noise mechanisms and to perform well on syntheticly generated videos. \
**Code/our_v2e/create_test_video.m** - Code that generates synthetic videos to be used by the v2e tool. 

## Threshold Estimation 
**Code/our_v2e/test_OffE_estimation.m** - this code containts the v2e tool and performing the OffE estimation algorithm suggested by Ziwei Wang et al at the paper "Event Camera Calibration of Per-Pixel Biased Contrast Threshold". \ 
We add analysis with our v2e tool to provide a ground truth for the threshold estimation results that can provide additional insight on the algorithm's performance in addition to the effect on video recostruction explored in the original paper.

## Threshold Mismatch Correction
**Code/threshold correction/main.m** - Implements a novel threshold mismatch correction algorithm based on sampling theory. The algorithm offers a new point of view on threshold mismatch correction, that of translating one threshold value to another rather than removing undesired events.\ 
This part is is stand alone and does not use the v2e and estimation code. 
