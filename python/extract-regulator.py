#!/usr/bin/python
import re

ode_file = "/home/spinicck/PhD/tmp/ode.txt"
csv_file = "/home/spinicck/PhD/Data/param-model/parameter/regulator.csv"
regex_gene = "dx\((\d+)\)"
regex_regulators = "\*\(S\((\d+),"
genes = ["Akt","AMPK","cMyc","HIF.1","mTOR","NOX","p53","PDK",
                "PI3K","PTEN","RAS","SOD","VEGF","GluT1","HK","G6PD.6PGD",
                "GPI","PFKFB2.3","PFK.1","ALD","TPI","GAPDH","PGK","PHGDH",
                "PGAM", "ENO","PKM2","PDH","ACC","LDH", "Glucose","G6P","F6P",
                "FBP", "G3P","DHAP","1,3BPG","3PG","2PG","PEP", "Pyruvate",
                "Lactate", "R5P","F2,6BP","Serine","Citrate","AMP","ADP", "ATP",
                "NAD","NADH","complex2","ROS"]

with open(ode_file, "r") as input:
    with open(csv_file, "w") as output:
        for line in input:
            gene = re.findall(regex_gene, line)
            all_reg_idx = re.findall(regex_regulators, line)
            csv_line = f'"{genes[int(gene[0])-1]}",'
            for i in all_reg_idx:
                csv_line = csv_line + f'"{ genes[int(i)-1] }",'
            csv_line = csv_line.rstrip(', ')
            csv_line = csv_line + "\n"
            output.write(csv_line)