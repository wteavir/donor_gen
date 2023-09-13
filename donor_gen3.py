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
                if number %2 != 0:
                    
            
                    insert_seq = feature.extract(insert_object)
                    #print(insert_seq)
                    
                    odd_backbone_object = SeqIO.read(odd_backbone, 'gb')
                    
                    #Finding the AscI & NotI location
                
                    NotI= "GCGGCCGC"
                    AscI= "GGCGCGCC"
                    Hx= "CGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCG"
                    NotI_position = odd_backbone_object.seq.find(NotI)
                    AscI_position = odd_backbone_object.seq.find(AscI)
                    
                    #removing Hx sequence 
                    Hx_position = odd_backbone_object.seq.find(Hx)
                    Hx_Start_End = [Hx_position,Hx_position+40]
                    
                    sequence_from_NotI = odd_backbone_object[NotI_position:]
                    AscI_NotI = [AscI_position + 8,NotI_position -1]
                    sequence_up_to_AscI = odd_backbone_object[:AscI_NotI[0]]
                    
                    
                    sequence_before_Hx = odd_backbone_object[:Hx_Start_End[0]]
                    donor_seq_features = sequence_before_Hx + AscI + insert_seq + sequence_from_NotI
                    
                    #Generating GBK file for combined sequence
                    gbk_file_name = f"{output_folder}/odd_donor_{label}.gbk"
                    donor_seq_features.annotations["molecule_type"] = "DNA"              
                    # Save the donor sequence as a GenBank file
                    SeqIO.write(donor_seq_features, gbk_file_name, "genbank")
                    
                    if cleaned_label == "fragment#1":
                        donor_seq_features = sequence_up_to_AscI + insert_seq + sequence_from_NotI
                        donor_seq_features.annotations["molecule_type"] = "DNA"  
                        SeqIO.write(donor_seq_features, gbk_file_name, "genbank")
                        
                    
                if number % 2 == 0:
                     
             
                     insert_seq = feature.extract(insert_object)
                     #print(insert_seq)
                     
                     even_backbone_object = SeqIO.read(even_backbone, 'gb')
                     
                     
             
                     #Finding the AscI & NotI location
                 
                     NotI= "GCGGCCGC"
                     AscI= "GGCGCGCC"
                     NotI_position = even_backbone_object.seq.find(NotI)
                     AscI_position = even_backbone_object.seq.find(AscI)
                     AscI_NotI = [AscI_position + 8,NotI_position -1]
                     
                     #Combining Sequence to AscI + Insert + Sequence from NotI 
                     sequence_up_to_AscI = even_backbone_object[:AscI_NotI[0]]
                     
                     sequence_from_NotI = even_backbone_object[NotI_position:]
                     donor_seq_features = sequence_up_to_AscI + insert_seq + sequence_from_NotI
                     donor_seq_features.annotations["molecule_type"] = "DNA"
                     
                     
                     
                     
                     #Generating GBK file for combined sequence
                     gbk_file_name = f"{output_folder}/even_donor_{label}.gbk"
                                   
                     # Save the donor sequence as a GenBank file
                     SeqIO.write(donor_seq_features, gbk_file_name, "genbank")
              
                                        
          
#donor_gen(backbone_file, odd_backbone,even_backbone, insert_file,output_folder)



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
