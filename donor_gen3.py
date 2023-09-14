#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 10:47:10 2023

@author: wteavir1
"""

from Bio import SeqIO
import os
import argparse
import re

#backbone_file = "/Users/wteavir1/Documents/Genbank_Files/Backbone"
#insert_file = "/Users/wteavir1/Documents/Genbank_Files/Backbone_Insert/mpapaya.gb"
#output_folder = "/Users/wteavir1/Documents/Output_folder/donor_gen"

odd_backbone = 'backbone_psl1139.gb'
even_backbone = 'psl1213-psl1194-kanr-apmr.gb'

NotI= "GCGGCCGC"
AscI= "GGCGCGCC"
Hx= "CGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCG"

def donor_gen (odd_backbone,even_backbone, insert_file,output_folder) :
    insert_object = SeqIO.read(insert_file, 'gb')
           
    #Isolating the feature with "Fragment" annotation
    features = insert_object.features
    for feature in features:
        qualifiers = feature.qualifiers
        if "label" in qualifiers: 
            label = qualifiers["label"][0]
            #Cleaning the Label 
            cleaned_label= label.replace(" ","").lower()
            match = re.search(r"fragment#?(\d+)",cleaned_label)
            if match:
                number = int(match.group(1))

                if number % 2 == 0:
                    backbone_object = SeqIO.read(even_backbone, 'gb')
                    gbk_file_name = f"{output_folder}/even_donor_{label}.gbk"

                else:
                    backbone_object = SeqIO.read(odd_backbone , 'gb')
                    gbk_file_name = f"{output_folder}/odd_donor_{label}.gbk"
                     
                insert_seq = feature.extract(insert_object)
                #print(insert_seq)
                
                #Finding the AscI & NotI location
                NotI_position = backbone_object.seq.find(NotI)
                AscI_position = backbone_object.seq.find(AscI)
                sequence_from_NotI = backbone_object[NotI_position:]
                sequence_up_to_AscI = backbone_object[:AscI_position + 8]

                if number % 2 != 0 and cleaned_label == "fragment#1":
                    #removing Hx sequence 
                    Hx_position = backbone_object.seq.find(Hx)
                    sequence_before_Hx = backbone_object[:Hx_position]

                    donor_seq_features = (sequence_before_Hx + AscI + 
                                insert_seq + sequence_from_NotI)

                donor_seq_features = sequence_up_to_AscI + insert_seq + sequence_from_NotI

                donor_seq_features.annotations["molecule_type"] = "DNA"  
                # Save the donor sequence as a GenBank file
                SeqIO.write(donor_seq_features, gbk_file_name, "genbank")
                        

if __name__ == '__main__':
    parser= argparse.ArgumentParser(description='Makes Donor sequence + features for every fragment')
    parser.add_argument('odd_backbone', help='file to the odd backbone sequence')
    parser.add_argument('even_backbone', help='file to the even backbone sequence')
    parser.add_argument('insert', help= 'file to the gene sequence')
    parser.add_argument('-o','--output', help='output folder directory for the donor sequences')
    args=parser.parse_args()
        
    donor_gen(args.odd_backbone,args.even_backbone, args.insert, args.output)
        
    if not os.path.isfile(args.odd_backbone):
        raise FileNotFoundError(f"GenBank file not found {args.odd_backbone}")
    if not os.path.exists(args.even_backbone):
        raise FileNotFoundError(f"Output Folder not found {args.even_backbone}")
    if not os.path.exists(args.insert):
        raise FileNotFoundError(f"Output Folder not found {args.insert}")
    if not os.path.exists(args.output):
        print(f"Creating output directory: {args.output}")
        os.makedirs(args.output)
    else:
        print(f"Output directory already exists: {args.output}")
