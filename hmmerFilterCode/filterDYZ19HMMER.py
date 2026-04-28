import pandas as pd
import numpy as np
import os
import collections
from collections import defaultdict
import gzip
import argparse as ap


def main():
    global args
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, help='What is the path to the input folder /my/input/folder/hmmerFiles/')
    parser.add_argument("-o", "--output", dest='output', required=True, help='Give the full path to output the dataframe')
    args = parser.parse_args()

    # File path to nhmmer output
    tblout_folder = str(args.input)
    
    goodFiles=['DYZ19_Yq.hmmer-table.txt.gz']
    # Dictionary to store the best hits per target coordinate
    
    header=['target_name','accession', 'query_name','accession2', 'hmmfrom', 'hmm_to', 'alifrom','ali_to','envfrom', 'env_to', 'sq_len', 'strand', 'E_value', 'score' ,'bias', 'description_of_target']
    
    # Parse the nhmmer tblout for each file and place into final blocks to build new dataframe that combines the 'good blocks'
    finalBlocks=[]
    for hmmerFile in os.listdir(tblout_folder):
        if '.hmmer-tblout.txt.gz' in hmmerFile:
            if "DYZ19" in hmmerFile:
                print(hmmerFile)
                best_hits = []
                with gzip.open(tblout_folder+"/"+hmmerFile, 'rt') as file:
                    for line in file:
                        if line.startswith("#"):
                            continue
                            
                        else:
                            best_hits.append(line.split())
                           
                file.close()
                
                tempDF = pd.DataFrame(data=best_hits, columns=header)
                tempDF['ali_to']=tempDF['ali_to'].astype(int)
                tempDF['alifrom']=tempDF['alifrom'].astype(int)
                tempDF['env_to']=tempDF['env_to'].astype(int)
                tempDF['envfrom']=tempDF['envfrom'].astype(int)
                tempDF['E_value']=tempDF['E_value'].astype(float)
                
                tempDF.sort_values(by=['target_name','alifrom'], inplace=True)
                
                for row in tempDF.index:
        
                    alifrom = tempDF.at[row,'alifrom']
                    alito = tempDF.at[row,'ali_to']
                    envfrom = tempDF.at[row,'envfrom']
                    envto = tempDF.at[row,'env_to']
                    
                    if int(envto) > int(envfrom):
                       continue
                    else:
                        tempDF.at[row,'env_to']=envfrom
                        tempDF.at[row,'envfrom']=envto
                        tempDF.at[row,'ali_to']=alifrom
                        tempDF.at[row,'alifrom']=alito
                
                #Measure the distance between hits
                tempDF['alignmentDistance']='TEMP'
                tempDF['envelopeDistance']='TEMP'
                flag=0
                for row in tempDF.index:
                    if flag==0:
                        tempDF.at[row,'alignmentDistance']=0
                        tempDF.at[row,'envelopeDistance']=0
                
                        previousali = int(tempDF.at[row,'ali_to'])
                        previousenv = int(tempDF.at[row,'env_to'])
                        flag+=1
                    else:
                        tempDF.at[row,'alignmentDistance']= int(tempDF.at[row,'alifrom'])-previousali
                        tempDF.at[row,'envelopeDistance']=int(tempDF.at[row,'envfrom'])-previousenv
                
                        previousali = int(tempDF.at[row,'ali_to'])
                        previousenv = int(tempDF.at[row,'env_to'])
        
        
                #Start putting together good blocks of hits
                goodBlocks={}
                blockNumber=0
                myflag2=0
                for row in tempDF.index:
                    if myflag2 ==0:
                        goodBlocks[blockNumber]={'rows':[row], 'Element':str(tempDF.at[row,'query_name'])}
                        myflag2+=1
                    else:
                        if abs(tempDF.at[row,'envelopeDistance'])<100 and (tempDF.at[row,'query_name'] == goodBlocks[blockNumber]['Element']):
                            goodBlocks[blockNumber]['rows'].append(row)
                        else:
                            blockNumber+=1
                            goodBlocks[blockNumber]={'rows':[row], 'Element':str(tempDF.at[row,'query_name'])}
        
                #Reform Dataframe
                for block in goodBlocks:
        
                    #Grab only those rows
                    tempDF2 = tempDF.loc[goodBlocks[block]['rows']].copy()
                    tempDF2.reset_index(inplace=True)
        
                    #make a total length dictionary to see what proportion of the block is zero evalue
                    totalLengthDict = {x:0 for x in range(min(tempDF2['envfrom']), max(tempDF2['env_to'])+1)}
        
                    for row in tempDF2.index:
                        if tempDF2.at[row,'E_value']==0.0:
                            for number in range(int(tempDF2.at[row,'envfrom']), int(tempDF2.at[row,'env_to'])+1):
                                totalLengthDict[number]=1
                        else:
                            continue
                        
                    totalSize = sum(totalLengthDict.values())/len(totalLengthDict)
        
                    #What is the standard deviation of the distances between hits
                    blockSTD = np.std(tempDF2.iloc[1:]['envelopeDistance'])
                
                    
                    myRows = [[rowvalue for rowvalue in tempDF2.loc[x]] for x in tempDF2.index]
                    contig = [x for x in tempDF2['target_name']][0]
                    element = [x for x in tempDF2['query_name']][0]
                    start = min([int(y) for y in tempDF2['envfrom']])
                    end = max([int(x) for x in tempDF2['env_to']])
                
                    if len(set([x for x in tempDF2['strand']])) >1:
                        orientation = 'Possible_Inversion'
                    else:
                        orientation = [x for x in tempDF2['strand']][0]
                        
                    finalBlocks.append([contig, element, orientation, start, end, len(totalLengthDict), len([numberone for numberone in totalLengthDict.values() if numberone !=0])/len(totalLengthDict), blockSTD, len(myRows), myRows])
                
                else:
                    continue
            else:
                continue
    
    #Build New Dataframe
    fullBlockDF = pd.DataFrame(data=finalBlocks).dropna()
    fullBlockDF.sort_values(by=[0,3,4], inplace=True)
    
    if len(fullBlockDF) ==1:
        finalDF = fullBlockDF[[0,1,2,3,4,8]].copy()
        finalDF.columns=['Contig','Element','Orientation','Start','End','Total_DYZ19_Copies']
        finalDF['ArrayBlocks_Combined']=1
        finalDF.reset_index().drop(columns=['index']).set_index('Contig', inplace=True)
        finalDF2 = finalDF.reset_index().drop(columns=['index']).set_index('Contig').copy()
        finalDF2.to_csv(str(args.output))
    else:
        #measure distance between blocks
        fullBlockDF['BlockDistance']='TEMP'
        flag=0
        for row in fullBlockDF.index:
            if flag==0:
                fullBlockDF.at[row,'BlockDistance']=0
        
                previousali = int(fullBlockDF.at[row,4])
                flag+=1
            else:
                fullBlockDF.at[row,'BlockDistance']= int(fullBlockDF.at[row,3])-previousali
        
                previousali = int(fullBlockDF.at[row,4])
    
    
        #make a final df allow some distance in case some weird MEI happens
        if len(set(fullBlockDF[0]))==1:
            fullBlockDF2 = fullBlockDF[fullBlockDF['BlockDistance']<=10000].copy()
            contig = list(fullBlockDF[0])[0]
            element = list(fullBlockDF[1])[0]
            if len(set(fullBlockDF[2]))==1:
                orientation = list(fullBlockDF[2])[0]
            else:
                orientation = 'Possible_Inversion_In_Array'
    
            start = min([int(x) for x in fullBlockDF[3]])
            end = max([int(x) for x in fullBlockDF[4]])
            totalElements = sum([int(x) for x in fullBlockDF[8]])
            
            finalDF = pd.DataFrame(data=[contig, element, orientation, start, end, totalElements, len(fullBlockDF)]).T
            finalDF.columns=['Contig','Element','Orientation','Start','End','Total_DYZ19_Copies', 'ArrayBlocks_Combined']
            finalDF.reset_index().drop(columns=['index']).set_index('Contig', inplace=True)
            finalDF2 = finalDF.reset_index().drop(columns=['index']).set_index('Contig').copy()
            finalDF2.to_csv(str(args.output))
    
        else:
            finalDFList=[]
            for contig in set(fullBlockDF[0]):
                fullBlockDF2 = fullBlockDF[(fullBlockDF['BlockDistance']<=10000) & (fullBlockDF[0]==contig)].copy()
                contig = list(fullBlockDF2[0])[0]
                element = list(fullBlockDF2[1])[0]
                if len(set(fullBlockDF2[2]))==1:
                    orientation = list(fullBlockDF2[2])[0]
                else:
                    orientation = 'Possible_Inversion_In_Array'
        
                start = min([int(x) for x in fullBlockDF2[3]])
                end = max([int(x) for x in fullBlockDF2[4]])
                totalElements = sum([int(x) for x in fullBlockDF2[8]])
                finalDFList.append([contig, element, orientation, start, end, totalElements, len(fullBlockDF2)])
                
            finalDF = pd.DataFrame(data=finalDFList, columns=['Contig','Element','Orientation','Start','End','Total_DYZ19_Copies', 'ArrayBlocks_Combined'])
            finalDF2 = finalDF.reset_index().drop(columns=['index']).set_index('Contig').copy()
            finalDF2.to_csv(str(args.output))

if __name__=="__main__":
    main()
