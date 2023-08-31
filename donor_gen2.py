#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 11:04:16 2023

@author: wteavir1
"""


import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import argparse
import os
from Bio import SeqIO
import re

#backbone_file = "/Users/wteavir1/Documents/Genbank_Files/Backbone"
#insert_file = "/Users/wteavir1/Documents/Genbank_Files/Backbone_Insert/mpapaya.gb"
#output_folder = "/Users/wteavir1/Documents/Output_folder/donor_gen"

#odd_backbone = 'backbone_psl1139.gb'
#even_backbone = 'psl1213-psl1194-kanr-apmr.gb'






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
    parser = argparse.ArgumentParser(description='Makes Donor sequence + features for every fragment')
    parser.add_argument('--odd_backbone', help='file to the odd backbone sequence')
    parser.add_argument('--even_backbone', help='file to the even backbone sequence')
    parser.add_argument('--insert', help='file to the gene sequence')
    parser.add_argument('--output', help='output folder directory for the donor sequences', default='output')
    args = parser.parse_args()

    if args.odd_backbone and args.even_backbone and args.insert and args.output:
        donor_gen(args.odd_backbone, args.even_backbone, args.insert, args.output)




def browse_file(entry):
    filename = filedialog.askopenfilename(filetypes=[("GenBank Files", "*.gb")])
    entry.delete(0, tk.END)
    entry.insert(0, filename)

def browse_output_folder(entry):
    foldername = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, foldername)
    

def execute_script():
    odd_backbone = entry_odd.get()
    even_backbone = entry_even.get()
    insert = entry_insert.get()
    output = entry_output.get()
    if odd_backbone and even_backbone and insert and output:
        donor_gen(odd_backbone, even_backbone, insert, output)
        messagebox.showinfo("Yayyyy", "You did it!")
    else:
        messagebox.showerror("Noooo", "Awwww! Please provide all input paths and output folder.")
    
    
root = tk.Tk()
root.title("Donor Sequence Generator for the Win")

# I did not write this interface button code, copied it from the internet sorry!
entry_odd = tk.Entry(root, width=40)
entry_odd.grid(row=0, column=1, padx=10, pady=5)
btn_browse_odd = tk.Button(root, text="Odd", command=lambda: browse_file(entry_odd))
btn_browse_odd.grid(row=0, column=2)

entry_even = tk.Entry(root, width=40)
entry_even.grid(row=1, column=1, padx=10, pady=5)
btn_browse_even = tk.Button(root, text="Even", command=lambda: browse_file(entry_even))
btn_browse_even.grid(row=1, column=2)

entry_insert = tk.Entry(root, width=40)
entry_insert.grid(row=2, column=1, padx=10, pady=5)
btn_browse_insert = tk.Button(root, text="Insert", command=lambda: browse_file(entry_insert))
btn_browse_insert.grid(row=2, column=2)

entry_output = tk.Entry(root, width=40)
entry_output.grid(row=3, column=1, padx=10, pady=5)
btn_browse_output = tk.Button(root, text="Output", command=lambda: browse_output_folder(entry_output))
btn_browse_output.grid(row=3, column=2)

# Call the execute_script function when all inputs are selected
entry_odd.bind("<FocusOut>", lambda event: execute_script())
entry_even.bind("<FocusOut>", lambda event: execute_script())
entry_insert.bind("<FocusOut>", lambda event: execute_script())
entry_output.bind("<FocusOut>", lambda event: execute_script())

btn_execute = tk.Button(root, text="Execute Script", command=execute_script)
btn_execute.grid(row=5, column=1, pady=10)

root.mainloop()


#Win Teavir

