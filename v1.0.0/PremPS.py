#!/usr/bin/python
# coding=utf-8
import os, sys, os.path, commands, re, getopt, datetime, time
from string import ascii_uppercase
from string import ascii_lowercase
from collections import defaultdict
import pandas as pd
import rpy2.robjects as robjects

ascii_cases = ascii_uppercase + ascii_lowercase + ''.join([str(i) for i in range(0,10)])  # if this case, the maximal inputchains should be 62.
r = robjects.r
r('''library(randomForest)''')

# set up path.
workdir = "/data/webservice/premps/"  # can change
pathpara = workdir + "inputfiles"
pathvmd = '/usr/local/bin/vmd'
pathmkdssp = '/usr/local/bin/mkdssp'
pathpsiblast = '/data/zhaofy/blast/psiblast'
pathblastdb = '/data/zhaofy/blastdb/nr'
pathrscript = '/usr/local/bin/Rscript'

# set up input and output file name for each job
jobid = ''
myopts, args = getopt.getopt(sys.argv[1:], "i:z")
for o, a in myopts:
    if o == '-i':
        jobid = a
jobpath = workdir + jobid
pathoutput = workdir + jobid + "out"  # can change

in_file = jobpath + '/' + jobid + '.input'
out_file = jobpath + '/' + jobid + '.sunddg'

os.system('mkdir %s' % pathoutput)

normal_format_pro = ['CYS','GLN','ILE','SER','VAL','MET','ASN','PRO','LYS','THR','PHE','ALA','HIS','GLY','ASP','LEU','ARG','TRP','GLU','TYR','PTR','TPO']
# map residue name three letters to one
map_three_one = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
                 "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
                 "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
                 "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
                 "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
                 "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "p",
                 "PTR": "y", "XLE": "X", "XAA": "X", "HSD": "H", "HID": "H",
                 "HSE": "H"}

map_one_three = {"G": "GLY", "A": "ALA", "S": "SER", "T": "THR", "C": "CYS",
                 "V": "VAL", "L": "LEU", "I": "ILE", "M": "MET", "P": "PRO",
                 "F": "PHE", "Y": "TYR", "W": "TRP", "D": "ASP", "E": "GLU",
                 "N": "ASN", "Q": "GLN", "H": "HIS", "K": "LYS", "R": "ARG",
                 "y": "PTR", "p": "TPO"}

map_surface = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,
               'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,
               'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,
               'Y':236.8,'V':164.5,'X':88.1}

# PON-tstab: Protein Variant Stability Predictor. Importance of Training Data Quality
hydrophobic = ['V','I','L','F','M','W','Y','C']
negatively_charged = ['D','E']
positively_charged = ['R','K']


# 从1A43.pdb中提取以ATOM开头的部分，加一列链的信息 A_1，写入1a43_p.pdb
def ProPDB1():
    pdball = []
    with open(in_file, 'r') as f:
        f.next()
        for line in f:
            ff = line.split('\t')
            pdb = ff[0].split('.')[0]  # 1AUD.pdb or 1AUD.pdb1
            suffix = ff[0].split('.')[1]
            usedchains = ff[1]
            if pdb not in pdball:
                pdball.append(pdb)
                ffpdb = open(pathoutput + '/' + pdb.lower() + '_p.pdb', 'w')
                try:
                    fpdb = open(jobpath + '/' + pdb.upper()+'.'+suffix, 'r')
                except:
                    fpdb = open(jobpath + '/' + pdb.lower()+'.'+suffix, 'r')
                if suffix != "pdb":
                    ST_play = False
                    for linepdb in fpdb:
                        if linepdb[0:5] == "MODEL":
                            CountModel = linepdb.split()[1]
                            ST_play = True
                            continue
                        if ST_play:
                            if linepdb[:4] == "ATOM" and linepdb[21] in usedchains and linepdb[17:20].strip() in normal_format_pro:
                                ffpdb.write("%s                  %s\n" % (
                                linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
                else: # .pdb
                    ST_play = True
                    for linepdb in fpdb:
                        line_list = re.split(r'\s+', linepdb)
                        if (line_list[0] == 'MODEL') and ST_play:
                            countmodel = line_list[1]
                            ST_play = False
                        if (line_list[0] == 'MODEL') and (line_list[1] != countmodel):
                            break
                        #if (line_list[0] == 'MODEL') and (line_list[1] == '2'):
                        #    break
                        if linepdb[:4] == 'ATOM' and linepdb[21] in usedchains and linepdb[17:20].strip() in normal_format_pro:
                            ffpdb.write("%s                  %s\n" % (linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
            else:
                continue
    ffpdb.close()
    fpdb.close()


# 删除1tm3_p.pdb中不在20个氨基酸内的残基，以及原子信息不同时包括['C','N','O','CA']的残基
# 首先生成1tm3_p_test.pdb, 再move成1tm3_p.pdb
def del_unknown_incomplete():
    residues = [k for k, v in map_three_one.iteritems() if v!= 'X']
    f = open(in_file, 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        pdbid = ff[0].split('.')[0].lower()
        pdbfile = open('{}/{}_p.pdb'.format(pathoutput,pdbid)).readlines()
        # delete unknow resdues
        pdbfile_del_unknown = [i for i in pdbfile if i[17:20] in residues]
        # delete incomplete residues
        final_row = pdbfile_del_unknown[-1]  # 最后一行
        last = ''          # 上一行的残基信息
        above = []         # 每一个残基的所有原子信息
        allresidues = []   # 所有残基的信息
        for row in pdbfile_del_unknown:
            if row[17:26] == last and row == final_row:            # when read final row，append it if equal to last
                above.append(row)
                atoms = [i[13:16].strip() for i in above if i[16]==' ' or i[16]=='A']          # ！！！16位可能是A，B表示状态， 不能取
                if set(['C','N','O','CA']).issubset(set(atoms)):
                    allresidues.append(above)
            elif row[17:26] == last and row != final_row:            # when read same residue, but not last row
                above.append(row)
            else:                            # when read different residue
                if len(above)>=4:
                    atoms = [i[13:16].strip() for i in above if i[16]==' ' or i[16]=='A']   # ！！！16位可能>是A，B表示状态， 不能取
                    if set(['C','N','O','CA']).issubset(set(atoms)):
                        allresidues.append(above)
                above = [row]
            last = row[17:26]
        # write out
        with open('{}/{}_p_test.pdb'.format(pathoutput,pdbid),'w') as fw:
            fw.write(''.join([y for x in allresidues for y in x]))
        break
    os.system('mv {}/{}_p_test.pdb {}/{}_p.pdb'.format(pathoutput,pdbid, pathoutput,pdbid))
    f.close()


# 从1a43out/1a43_p.pdb中抓取每条链的所有信息，从而为野生型蛋白质的每条链产生一个pdb文件， 1A43_A_1.pdb, 1A43_B_1.pdb.....
def splitchain():
    pdball = []
    with open(in_file, 'r') as f:
        f.next()
        for line in f:
            ff = line.split('\t')
            pdbfile = ff[0]  # 1A43.pdb
            pdb = pdbfile.split('.')[0]  # 1A43
            partner1 = ff[1].split('.')
            if pdbfile not in pdball:
                for chain in list(partner1):
                    f1 = open(pathoutput + '/' + pdb.lower() + '_p.pdb', 'r')
                    fw = open(pathoutput + '/' + pdb + '_' + chain + '.pdb', 'w')
                    for line1 in f1:
                        if chain == line1[72:].strip():
                            fw.write(line1)
                pdball.append(pdbfile)


# processing pdb files and prepair command files
# produce 1A43_A_1.seq, 1A43_CH1.pdb and 1A43_CH2.pdb by 1A43_A_1.pdb and 1A43_B_1.pdb.
# write files：1a43.input.cleaned
def CleanPdb():
    first_line = file(in_file, 'r').readlines()[0][:-1]
    fw = open(in_file + ".cleaned", "w")  # 1a43.input.cleaned
    fw.write("%s\t%s\t%s\t%s\n" % (first_line, "PDBid", "NewChains", "Mutation_cleaned"))

    second_line = file(in_file, 'r').readlines()[1]  # file函数类似于open函数
    ff = second_line.split("\t")
    pdb = ff[0].split(".")[0]  # 1A43
    partner1 = ff[1].split(".")  # [A_1]

    # map chain, renumber chain
    mapchainarray = []
    counti = 0
    for chains in list(partner1):
        cc = (chains, ascii_cases[counti])
        mapchainarray.append(cc)  # mapchainarray = [('A_1', 'A'), ('B_1', 'B'), ('C_1', 'C')]
        counti += 1
        mapchaindict = dict(iter(mapchainarray))  # mapchaindict = {'A_1': 'A', 'B_1': 'B', 'C_1':C}

    newpartner1 = ''
    for chains in list(partner1):
        newpartner1 += mapchaindict[chains]

    countchain = 0
    mut_clean = ""
    for chains in list(partner1):
        fvar = open(pathoutput + "/" + pdb + "_" + chains + ".var", "w")  # 1a43out/1A43_A_1.var
        countchain += 1
        count = 1
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")  # 1a43out/1A43_A_1.pdb
        line1 = fpdb.readlines()[0]  # 第一行内容
        resname = line1[17:20].strip()  # ALA  .or.  A/T/C/G  .or.  DA/AU/DG/AC
        resnum = line1[22:27].strip()  # 1
        # 写入1A43_A_1.seq, 1A43_A_1.var的内容
        mutchainall = []
        f = open(in_file, 'r')
        f.next()
        for line in f:
            ff = line.split("\t")
            mut = ff[3].strip()  # G9A
            mutchain = ff[2]  # A_1
            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                fseq = open(pathoutput + "/" + pdb + "_" + mutchain + ".seq", "w")  # 1a43out/1A43_A_1.seq
                fseq.write("%s %s\n" % (">", pdb + mutchain))  # > 1A43_A_1
                fseq.write("%s" % (map_three_one[resname]))  # map_three_one[ALA]
            mutchainall.append(mutchain)
            if str(chains) == str(mutchain):
                if str(resnum) == str(mut[1:-1]) and str(map_three_one[resname]) == str(mut[0:1]):
                    mut_clean = str(map_three_one[resname]) + mapchaindict[chains] + str(count) + mut[-1:]
                    fw.write("%s\t%s\t%s\t%s\n" % (line[:-1], pdb, newpartner1, mut_clean))
                    fvar.write("%s\n" % (str(map_three_one[resname]) + str(count) + mut[-1:]))
        f.close()

        fwpdb = open(pathoutput + "/" + pdb + "_CH" + str(countchain) + ".pdb", "w")
        fwpdb.write(
            "%s %s%s   %s %s" % (line1[0:16], line1[17:21], mapchaindict[chains], str(count), line1[27:]))
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        for linepdb in fpdb.readlines()[1:]:
            if linepdb[16:17] == " " or linepdb[16:17] == "A":
                resnamepdb = linepdb[17:20].strip()
                resnumpdb = linepdb[22:27].strip()
                if resnamepdb == resname and resnumpdb == resnum:
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 10 <= count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 100 <= count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 1000 <= count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                else:
                    count += 1
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 10 <= count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 100 <= count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if 1000 <= count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (
                            linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))

                    mutchainall = []
                    f = open(in_file, 'r')
                    f.next()
                    for line in f:
                        ff = line.split("\t")
                        mut = ff[3]  # G9A
                        mutchain = ff[2]  # A_1
                        if str(chains) == str(mutchain) and mutchain not in mutchainall:
                            fseq.write("%s" % (map_three_one[resnamepdb]))
                        mutchainall.append(mutchain)
                        if str(resnumpdb) == str(mut[1:-1]) and str(map_three_one[resnamepdb]) == str(mut[0:1]) and str(
                                chains) == str(mutchain):
                            mut_clean = str(map_three_one[resnamepdb]) + mapchaindict[chains] + str(count) + mut[-1:]
                            fw.write(
                                "%s\t%s\t%s\t%s\n" % (line[:-1], pdb, newpartner1, mut_clean))
                            fvar.write("%s\n" % (str(map_three_one[resnamepdb]) + str(count) + mut[-1:]))

                    f.close()

                resname = linepdb[17:20].strip()
                resnum = linepdb[22:27].strip()

        fpdb.close()
        fwpdb.close()

    fseq.close()
    fvar.close()
    fw.close()


# produce renumbered and processed 1A43.pdb wild-type pdb structure
def wtpdb():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[6]  # 1A43
            if pdb not in pdball:
                commands.getoutput('cat %s/%s_CH*.pdb > %s/%s.pdb' % (pathoutput, pdb, pathoutput, pdb))
                pdball.append(pdb)


 # calculate secondary structure  with DSSP using wild type crystal structure of mutchain, 1A43_A_1.pdb, produce 1A43.dssp
def dssp():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.strip().split('\t')
            pdb = ff[6]  # 1A43
            if pdb not in pdball:
                pdball.append(pdb)
                os.system('%s %s/%s.pdb %s/%s.dssp' % (pathmkdssp,pathoutput, pdb, pathoutput, pdb))


# hydrophobicity_scales.txt
def hydrophobicity_scales():
    hydrophobicity_scales = pd.read_csv(pathpara + '/hydrophobicity_scales.txt', header=0, index_col=0,sep = '\t')
    with open('{}/{}_hydroscales.txt'.format(pathoutput, jobid), 'w') as fw, open(in_file + '.cleaned') as f:
        f.next()
        fw.write('{}\n'.format('\t'.join(['PDB_ID', 'Mutation_PDB', 'dhydroscales_omh'])))
        for line in f:
            ff = line.strip().split('\t')
            pdb = ff[6]  # 1A43
            Mutation_cleaned = ff[-1]  # G9A
            wild_residue = Mutation_cleaned[0]  # R
            mut_residue = Mutation_cleaned[-1]  # A
            hydroscales_omh_wt = str(hydrophobicity_scales.loc[map_one_three[wild_residue], 'OMH'])
            hydroscales_omh_mut = str(hydrophobicity_scales.loc[map_one_three[mut_residue], 'OMH'])
            dhydroscales_omh = str(float(hydroscales_omh_mut) - float(hydroscales_omh_wt))
            fw.write('{}\n'.format('\t'.join([pdb, Mutation_cleaned, dhydroscales_omh])))


def run_provean():
    mutchainall = []
    with open(in_file + ".cleaned", 'r') as f1:
        f1.next()
        for line in f1:
            ff = line.strip().split('\t')
            mutchain = ff[2]
            pdb = ff[6]
            mutation_cleaned = ff[-1]
            if mutchain not in mutchainall:
                mutchainall.append(mutchain)
                os.system('provean.sh -q %s/%s_%s.seq -v %s/%s_%s.var > %s/provean_%s_%s.out --num_threads 30' % (pathoutput,pdb,mutchain,pathoutput,pdb,mutchain,pathoutput,pdb,mutchain))


# PON-stab
def neighbor_res():
    with open(pathoutput+'/'+jobid+'.neighbor','w') as fw:
        fw.write('\t'.join(['PDBfile','Chains','MutChain','Mutation_PDB','Result_Id','isPI','PDBid','NewChains','Mutation_cleaned',
            'hydrophobic_wt_2d','charged_wt_2d'])+'\n')
        with open(in_file+'.cleaned') as f:
            f.next()
            for line in f:
                pdb = line.strip().split('\t')[6].lower()
                mutchain = line.strip().split('\t')[2]
                mut = line.strip().split('\t')[-1].lower()
                ### 2d
                seq_wt = open(pathoutput+'/'+pdb.upper()+'_'+mutchain+'.seq').readlines()[1]
                mutloc = int(mut[2:-1])
                if mutloc-11 >=0:
                    istart = mutloc-11-1
                else:
                    istart = 0
                if mutloc+11 <= len(seq_wt)-1:
                    iend = mutloc+11
                else:
                    iend = len(seq_wt)
                # wt
                hydrophobic_wt,negatively_charged_wt,positively_charged_wt = 0,0,0
                for i in range(istart,iend):
                    if seq_wt[i] in hydrophobic:
                        hydrophobic_wt += 1
                    if seq_wt[i] in negatively_charged:
                        negatively_charged_wt += 1
                    if seq_wt[i] in positively_charged:
                        positively_charged_wt += 1
                charged_wt_2d = negatively_charged_wt+positively_charged_wt

                fw.write(line.strip()+'\t'+'\t'.join([str(hydrophobic_wt),str(charged_wt_2d)])+'\n')


# SOLart
def solart():
    # 20 AA all classify(20*3)
    allclass = []
    for i in map_surface.keys():
        allclass = allclass+['buried_'+i,'mod_buried_'+i,'exposed_'+i]

    with open(in_file+'.cleaned') as f, open(pathoutput+'/'+jobid+'.solart','w') as fw:
        f.next()
        fw.write('PDBid\t%s\n' %('\t'.join(['buried_L','buried_charged_sum','buried_aromatic'])))
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[-3]
            chains = linelist[-2]
            daac = defaultdict(list)

            # Denominator: len_protein
            count_pro = set()
            with open(pathoutput + '/' + pdb + '.pdb', 'r') as f:
                for row in f:
                    count_pro.add((row[72:-1], row[22:27].strip()))
                pro_len = len(count_pro)

            # wt
            for chain in chains:
                with open(pathoutput+'/'+pdb+'.dssp') as fdssp:
                    for ldssp in fdssp:
                        if re.match(r'^\s+\d+\s+\d+ {}'.format(chain), ldssp):
                            res = ldssp[13]
                            resnum = ldssp[5:10].strip()
                            reschain = ldssp[11]

                            resacc = ldssp[35:39].strip()
                            resacc_tri = float(resacc)/map_surface[res]
                            if resacc_tri <= 0.2:
                                resloc = 'buried'
                            elif resacc_tri < 0.5:
                                resloc = 'mod_buried'
                            else:
                                resloc = 'exposed'

                            daac[resloc+'_'+res].append(reschain+'_'+resnum)


            # Amino acid composition
            # 20中aa组分
            aacdic = defaultdict(float)
            for aac in allclass:
                if aac in daac.keys():
                    aacdic[aac] = float(len(daac[aac]))/pro_len
                else:
                    aacdic[aac] = 0 
            # classify
            aacresult = [aacdic['buried_L'],
                         aacdic['buried_K']+aacdic['buried_R']+aacdic['buried_D']+aacdic['buried_E'],  # charged sum
                         aacdic['buried_F']+aacdic['buried_W']+aacdic['buried_Y']]  # aromatic
            aacresult = [str(i ) for i in aacresult]
            fw.write(pdb+'\t'+'\t'.join(aacresult)+'\n')
            break


# 位置特异性打分矩阵
def run_PSSM():
    with open(in_file+'.cleaned') as f:
        f.next()
        pdball = []
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[6]
            mutchain = linelist[2]
            inputseq_wt = pathoutput+'/'+pdb+'_'+mutchain+'.seq' # !!!
            output_wt = pathoutput+'/'+pdb+'_'+mutchain
            # run PSSM
            if pdb+mutchain not in pdball:
                pdball.append(pdb+mutchain)
                os.system('{} -query {} -db {} -num_iterations 3 -out_ascii_pssm {}.pssm'.format(pathpsiblast,inputseq_wt,pathblastdb,output_wt))


# get result form PSSM
def get_PSSM():
    dindex = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
    with open(in_file+'.cleaned') as f,open(pathoutput+'/'+jobid+'.pssm','w') as fw:
        title = f.next()
        fw.write(title.strip()+'\tposition_score_wtmut\n') 
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[6]
            mutchain = linelist[2]
            mut = linelist[-1]
            # wild
            with open(pathoutput+'/'+pdb+'_'+mutchain+'.pssm') as fpssm:
                for ffp in fpssm:
                    if re.match(r'\s+'+mut[2:-1]+'\s+'+mut[0], ffp):  # !!!
                        position_score_wtmut = ffp.strip().split()[dindex[mut[-1]]+2]

            # write
            fw.write(line.strip()+'\t'+position_score_wtmut+'\n')


# get all features to 1a43.input.cleaned.outdata
def getenergy():
    fw = open(in_file + ".cleaned.outdata", 'w')
    first_line = file(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            first_line, "location", "GXG_wt", "ACC_wt", 'provean', 'dhydroscales_omh', 'hydrophobic_wt_2d', 'charged_wt_2d',
            'buried_L', 'buried_charged_sum', 'buried_aromatic', 'position_score_wtmut'))

    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[6].lower()  # 1a43
        mutchain_nocleaned = ff[2]  # A_1
        mut = ff[8].lower()  # ga9a
        mutchain = mut[1:2]  # a
        resnum = mut[2:-1]  # 9
        protname = pdb + '_' + mut  # 1a43_ga9a
        pdball = []

        gxg_wt = map_surface[mut[0].upper()]

        # location
        with open(pathoutput + '/' + pdb.upper()+'.dssp', 'r') as fdssp:
            for ldssp in fdssp.readlines():
                if re.match(r'^\s+\d+\s+\d+ {}'.format(mutchain.upper()), ldssp):
                    if ldssp[5:10].strip() == mut[2:-1] :
                        ACC_wt = ldssp[35:39].strip()
        if float(ACC_wt)/gxg_wt <= 0.2:
            location = 'COR'
        else:
            location = 'SUR'

        # add hydro_scales of wild and mutant residue
        with open(pathoutput + '/{}_hydroscales.txt'.format(jobid)) as fhydro:
            fhydro.next()
            for lhydro in fhydro:
                lhydro_list = lhydro.strip().split('\t')
                if mut.upper() == lhydro_list[1]:
                    hydroscales = '\t'.join(lhydro_list[2:])

        # add pon-stab
        with open(pathoutput + '/{}.neighbor'.format(jobid)) as fnei:
            fnei.next()
            for lnei in fnei:
                lnei_list = lnei.strip().split('\t')
                if mut.upper() == lnei_list[8]:
                    neighbor_res = '\t'.join(lnei_list[9:])

        # solart
        with open(pathoutput + '/{}.solart'.format(jobid)) as fsol:
            fsol.next()
            for lsol in fsol:
                lsol_list = lsol.strip().split('\t')
                sol_ratio = '\t'.join(lsol_list[1:])

        # PSSM
        with open(pathoutput + '/{}.pssm'.format(jobid)) as fpssm:
            fpssm.next()
            for lpssm in fpssm:
                lpssm_list = lpssm.strip().split('\t')
                if mut.upper() == lpssm_list[8]:
                    pssm_out = '\t'.join(lpssm_list[9:])

        # add provean score
        ST_play = False
        with open(pathoutput+'/provean_'+pdb.upper()+'_'+mutchain_nocleaned+'.out','r') as fprovean:
            for ffprovean in fprovean:
                if ffprovean[2:11] == "VARIATION":
                    ST_play = True
                    continue
                if ST_play:
                    ffp = ffprovean.strip().split("\t")
                    if ffp[0] == (mut[0]+mut[2:]).upper():
                        provean_score = ffp[1]

        # write final file
        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                line[:-1], location, str(gxg_wt), str(ACC_wt), provean_score, hydroscales, neighbor_res, sol_ratio, pssm_out)) 

    f.close()
    fw.close()


# using our fitting models.
def Prediction():
    outdata = in_file + '.cleaned.outdata'
    robjects.globalenv["outdata"] = outdata
    robjects.globalenv["workdir"] = workdir
    r('''test = read.table(outdata,header=T,sep="\t")''')
    r('''filename = paste(workdir, 'inputfiles/PremPS_new.RData',sep = '')''')
    r('''load(file = filename)''')
    PredR = r('''predict(premps.rf,test)''')
    cutoff = 1 
    robjects.globalenv["PredR"] = PredR
    predDDG = r('''PredR''')

    first_line = file(in_file, 'r').readlines()[0][:-1]
    fw = open(out_file, "w")
    fw.write("%s\t%s\t%s\t%s\t%s\n" % (first_line, "SunDDG", "Location", "Deleterious", "Confidence"))

    f = open(outdata, 'r')
    _unused = f.next()
    count = 0
    for line in f:
        ff = line.split("\t")
        if abs(float(predDDG[count])) >= cutoff:
            Deleterious = "yes"
            Confidence = "Low"
        else:
            Deleterious = "no"
            Confidence = "Low"
        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%3.2f\t%s\t%s\t%s\n" % (
            ff[0], ff[1], ff[2], ff[3], ff[4], ff[5], predDDG[count], ff[9], Deleterious, Confidence))
        count += 1
    f.close()
    fw.close()


# Add RF contribution
def rf_contribution():
    outdata = in_file + '.cleaned.outdata'
    template = open(pathpara+'/PremPS_forestFloor_new.R').read()
    result = template.replace('test_outfeature',outdata).replace('test_outcontribution',outdata+'.contribution').replace('test_sunddg',out_file)
    with open(pathoutput+'/'+jobid+'_forestFloor.R','w') as fw:
        fw.write(result)
    os.system('%s %s' % (pathrscript,pathoutput+'/'+jobid+'_forestFloor.R'))


# produce psf and pdb files of wild-type with vmd.
# produce 1A43_vmd.psf and 1A43_vmd.pdb using 1A43_CH1.pdb and 1A43_CH2.pdb.
def vmd_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[6]
        partner1 = ff[7]
        NumChain = int(len(partner1))
        if pdb not in pdball:
            vmd_pdb = template.replace('protname', pdb).replace('NumChain', str(NumChain)).replace('pathinput', pathpara).replace('pathoutput', pathoutput)
            with open(pathoutput + '/vmd_' + pdb + '.pgn', 'w') as fw:
                fw.write(vmd_pdb)
            commands.getoutput('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd,pathoutput, pdb))
            pdball.append(pdb)
            # change HSD to HIS
            os.system('sed -e "s/HSD/HIS/g" %s/%s_vmd.pdb > %s/%s_vmd_temp.pdb' % (pathoutput, pdb, pathoutput, pdb))
            os.system('mv %s %s' % (pathoutput+'/'+pdb+'_vmd_temp.pdb',pathoutput+'/'+pdb+'_vmd.pdb'))
        else:
            continue
    f.close()


# make inputfile for foldx4 with mutchain for calculating folding free energy
# produce individual_list_1A43_A_GA9A.txt,  foldx_buildmodel_1A43_A_GA9A.txt,   1A43_A_Repair_GA9A.pdb
def inputfoldx():
    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        mut = ff[8][:-1]  # GA9A
        pdb = ff[6]  # 1A43
        with open('individual_list_' + pdb + '_' + mut + '.txt', 'w') as fic:
            fic.write('%s' % (mut + ';'))
        with open('foldx_buildmodel_' + pdb + '_' + mut + '.txt', 'w') as fsc:
            fsc.write('command=BuildModel\npdb=%s\nmutant-file=%s' % (
                pdb + '_' + mut + '.pdb', 'individual_list_' + pdb + '_' + mut + '.txt'))
        commands.getoutput("cp %s/%s.pdb %s_%s.pdb" % (pathoutput, pdb, pdb, mut))
    f.close()


# build model with foldx, produce Dif_1A43_A_Repair_GA9A.fxout
def runfoldx_mut():
    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        mut = ff[8][:-1]  # GA9A
        pdb = ff[6]  # 1A43

        # run buildmodel for pdbchain, partner1, complex
        os.system('./foldx -f foldx_buildmodel_%s_%s.txt' % (pdb, mut))

        # process outputfile produce by complex
        commands.getoutput("rm WT_%s_%s_1.pdb" % (pdb, mut))
        commands.getoutput("rm individual_list_%s_%s.txt" % (pdb, mut))
        commands.getoutput("rm foldx_buildmodel_%s_%s.txt" % (pdb, mut))
        commands.getoutput("rm Average_%s_%s.fxout" % (pdb, mut))
        commands.getoutput("rm Raw_%s_%s.fxout" % (pdb, mut))
        commands.getoutput("rm PdbList_%s_%s.fxout" % (pdb, mut))
        commands.getoutput("rm %s_%s.fxout" % (pdb, mut))
        commands.getoutput("mv Dif_%s_%s.fxout %s/" % (pdb, mut,pathoutput))
        commands.getoutput("rm %s_%s.pdb" % (pdb, mut))
        commands.getoutput("mv %s_%s_1.pdb %s/%s_%s.pdb" % (pdb, mut, pathoutput, pdb, mut))

    f.close()


# split chains and produce pdb files for each chain of mutant pdb, which are used for VMD.
# produce 1A43_G9A_CH1.pdb and 1A43_G9A_CH2.pdb by 1A43_G9A.pdb
def splitchain_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[6]
        partner1 = ff[7]
        mut = ff[8][:-1]
        pdb = pdb + '_' + mut
        countchain = 1
        for chains in list(partner1):
            os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s.pdb' % (
                chains, pathoutput, pdb, pathoutput, pdb, 'CH' + str(countchain)))
            countchain += 1
    f.close()


# produce psf and pdb files of mutant with vmd 
def vmd_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[6]
        partner1 = ff[7]
        mut = ff[8][:-1]
        result_id = ff[4]
        protname = pdb + '_' + mut
        NumChain = int(len(partner1))
        vmd_pdb = template.replace('protname', protname).replace('NumChain', str(NumChain)).replace('pathinput', pathpara).replace('pathoutput', pathoutput)
        with open(pathoutput + '/vmd_' + protname + '.pgn', 'w') as fw:
            fw.write(vmd_pdb)
        commands.getoutput('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, protname))
        with open(pathoutput+'/test.txt','w') as ftitle:
            ftitle.write("REMARK THIS MUTANT(" + pdb + '-' + ff[2] + '-' + ff[3] + ") STRUCTURE IS PRODUCED BY FOLDX\n")
            ftitle.write("REMARK DATE:" + "  " + time.strftime("%x") + "  " + time.strftime("%X") + "   CREATED BY SERVER: PREMPS\n")
            ftitle.write("REMARK REFERENCE: PLEASE CITE *** \n")
        os.system('cat %s %s > %s' % (pathoutput+'/test.txt', pathoutput+'/'+protname+'_vmd.pdb',pathoutput+'/'+pdb+'_'+result_id+'_vmd.pdb'))
        os.system('rm %s' % (pathoutput+'/test.txt'))
        os.system('sed -e "s/HSD/HIS/g" %s/%s_vmd.pdb > %s/%s_vmd_temp.pdb' % (pathoutput, pdb+'_'+result_id, pathoutput, pdb+'_'+result_id))
        os.system('mv %s %s' % (pathoutput+'/'+pdb+'_'+result_id+'_vmd_temp.pdb',pathoutput+'/'+pdb+'_'+result_id+'_vmd.pdb'))
    f.close()


def main():
    ProPDB1()
    del_unknown_incomplete()
    splitchain()
    CleanPdb()
    wtpdb()
    dssp()
    hydrophobicity_scales()
    run_provean()
    neighbor_res()
    solart()
    run_PSSM()
    get_PSSM()
    getenergy()
    Prediction()
    rf_contribution()
    vmd_wt()
    # run mutant
    ispi = file(in_file, 'r').readlines()[1].strip()[-1]
    if ispi == '1':
        inputfoldx()
        runfoldx_mut()
        splitchain_mut()
        vmd_mut()


#main
if __name__ == '__main__':
    main()
