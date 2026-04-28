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
    
    goodFiles=['TSPY.hmmer-table.txt.gz']
    # Dictionary to store the best hits per target coordinate
    
    header=['target_name','accession', 'query_name','accession2', 'hmmfrom', 'hmm_to', 'alifrom','ali_to','envfrom', 'env_to', 'sq_len', 'strand', 'E_value', 'score' ,'bias', 'description_of_target']
    
    # Parse the nhmmer tblout for each file and place into final blocks to build new dataframe that combines the 'good blocks'
    finalBlocks=[]
    for hmmerFile in os.listdir(tblout_folder):
        if '.hmmer-tblout.txt.gz' in hmmerFile:
            if "TSPY" in hmmerFile:
                print(hmmerFile)
                best_hits = []
                with gzip.open(tblout_folder+'/'+hmmerFile, 'rt') as file:
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
                        if abs(tempDF.at[row,'envelopeDistance'])<100 and tempDF.at[row,'query_name'] == goodBlocks[blockNumber]['Element']:
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
    finalDF = pd.DataFrame(data=finalBlocks)
    finalDF.sort_values(by=[0,3,4], inplace=True)
    
    goodTSPYRows=[]
    for row in finalDF.index:
        if finalDF.at[row,5] <10000:
            continue
        else:
            goodTSPYRows.append(row)
    goodTSPY1 = finalDF.loc[goodTSPYRows].copy()
    reCheckTSPY = finalDF.loc[[x for x in finalDF.index if x not in goodTSPYRows]].copy()

    individualTSPYHits=[]
    for grouping in goodTSPY1[9]:
        for tspyHit in grouping:
            if float(tspyHit[14])<1000:
                continue
            else:
                individualTSPYHits.append(tspyHit)
    allGoodTSPYHitDF = pd.DataFrame(data=individualTSPYHits).drop(columns=[2,4,7,8,11,16])
    
    refinedHits=[]
    mergedRows=[]
    for row in allGoodTSPYHitDF.index:
    
        #If we already merged this row just ignore it
        if row in mergedRows:
            continue
            
        else:
    
            #If this is the last row and it hasnt been merged then just add it
            if int(row) == len(allGoodTSPYHitDF)-1:
                
                refinedHits.append([x for x in allGoodTSPYHitDF.loc[row]])
    
    
            #This is the real test
            else:
                
                #If this row and the next are on the same contig then test it
                if allGoodTSPYHitDF.at[row,1] == allGoodTSPYHitDF.at[row+1,1]:
        
                    #If this row and the next are not overlapping at all just add this row
                    if int(allGoodTSPYHitDF.at[row,18]) >0 and int(allGoodTSPYHitDF.at[row+1, 18])>0:
                        
                        refinedHits.append([x for x in allGoodTSPYHitDF.loc[row]])
    
                    #If there is some overlap in the annotation
                    else:
            
                        if int(allGoodTSPYHitDF.at[row,5]) <50 and int(allGoodTSPYHitDF.at[row+1,6]) >20000 and (sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][14]])<25000):
                            mergedRows.append(row)
                            mergedRows.append(row+1)
            
                            refinedHits.append(['noNumber', allGoodTSPYHitDF.at[row,1], 'TSPY', min([int(b) for b in allGoodTSPYHitDF.loc[[row, row+1]][5]]),max([int(a) for a in allGoodTSPYHitDF.loc[[row, row+1]][6]]),
                                                min([int(y) for y in allGoodTSPYHitDF.loc[[row, row+1]][9]]),max([int(z) for z in allGoodTSPYHitDF.loc[[row, row+1]][10]]),
                                                allGoodTSPYHitDF.at[row,12], sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][13]]), sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][14]]), 'MergedRow', allGoodTSPYHitDF.at[row,17],
                                                allGoodTSPYHitDF.at[row,18]])
    
                        elif int(allGoodTSPYHitDF.at[row,6]) > 20000 and int(allGoodTSPYHitDF.at[row+1,5]) <50 and (sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][14]])<25000):
                            mergedRows.append(row)
                            mergedRows.append(row+1)
            
                            refinedHits.append(['noNumber', allGoodTSPYHitDF.at[row,1], 'TSPY', min([int(b) for b in allGoodTSPYHitDF.loc[[row, row+1]][5]]),max([int(a) for a in allGoodTSPYHitDF.loc[[row, row+1]][6]]),
                                                min([int(y) for y in allGoodTSPYHitDF.loc[[row, row+1]][9]]),max([int(z) for z in allGoodTSPYHitDF.loc[[row, row+1]][10]]),
                                                allGoodTSPYHitDF.at[row,12], sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][13]]), sum([float(d) for d in allGoodTSPYHitDF.loc[[row, row+1]][14]]), 'MergedRow', allGoodTSPYHitDF.at[row,17],
                                                allGoodTSPYHitDF.at[row,18]])
                        
                        else:
                            refinedHits.append([x for x in allGoodTSPYHitDF.loc[row]])
    
    
    
                
                            
                #if the next row isnt on the same contig then just add this row
                else:
                    refinedHits.append([x for x in allGoodTSPYHitDF.loc[row]])
            
    refinedDF = pd.DataFrame(data=refinedHits).drop(columns=[0])
    newGroups=[]
    tempGroup=[]
    flag=0
    for row in refinedDF.index:
        if flag==0:
            tempGroup.append(row)
            ending = int(refinedDF.at[row,6])
            contig = str(refinedDF.at[row,1])
            flag+=1
        else:
            if int(refinedDF.at[row,5])-ending < 10000 and contig == str(refinedDF.at[row,1]):
                tempGroup.append(row)
                ending = int(refinedDF.at[row,6])
                contig = str(refinedDF.at[row,1])
            else:
                newGroups.append(tempGroup)
                tempGroup=[]
                tempGroup.append(row)
                ending = int(refinedDF.at[row,6])
                contig = str(refinedDF.at[row,1])
                
    if len(tempGroup)>0:
        newGroups.append(tempGroup)
    
    # Find the shortest list
    shortest_list = min(newGroups, key=len)
    if len(shortest_list)==2:
        for row in shortest_list:
            if abs(int(refinedDF.at[row,4])-int(refinedDF.at[row,3]))>18000:
                refinedDF.at[row,2]='TSPY2'
            else:
                continue
    else:
        pass
    
    numberList = []
    for leftOverTSPY in reCheckTSPY[9]:
        for hit in leftOverTSPY:
            numberList.append(hit[0])
    redoDF = tempDF.loc[numberList].copy()
    #Measure the distance between hits
    redoDF['alignmentDistance']='TEMP'
    redoDF['envelopeDistance']='TEMP'
    flag=0
    for row in redoDF.index:
        if flag==0:
            redoDF.at[row,'alignmentDistance']=0
            redoDF.at[row,'envelopeDistance']=0
    
            previousali = int(redoDF.at[row,'ali_to'])
            previousenv = int(redoDF.at[row,'env_to'])
            flag+=1
        else:
            redoDF.at[row,'alignmentDistance']= int(redoDF.at[row,'alifrom'])-previousali
            redoDF.at[row,'envelopeDistance']=int(redoDF.at[row,'envfrom'])-previousenv
    
            previousali = int(redoDF.at[row,'ali_to'])
            previousenv = int(redoDF.at[row,'env_to'])
    
    #Start putting together good blocks of hits
    goodBlocks2={}
    blockNumber=0
    myflag2=0
    for row in redoDF.index:
        if myflag2 ==0:
            goodBlocks2[blockNumber]={'rows':[row], 'Element':str(redoDF.at[row,'query_name']) ,'hmmMatchTotal': int(redoDF.at[row,'env_to'])-int(redoDF.at[row,'envfrom'])}
            myflag2+=1
        else:
            if abs(redoDF.at[row,'envelopeDistance'])<750 and (redoDF.at[row,'query_name'] == goodBlocks2[blockNumber]['Element']):
                goodBlocks2[blockNumber]['rows'].append(row)
                addNumberHMM = int(redoDF.at[row,'hmm_to'])-int(redoDF.at[row,'hmmfrom'])
                goodBlocks2[blockNumber]['hmmMatchTotal']+=addNumberHMM
            else:
                blockNumber+=1
                goodBlocks2[blockNumber]={'rows':[row], 'Element':str(redoDF.at[row,'query_name']) ,'hmmMatchTotal': int(redoDF.at[row,'env_to'])-int(redoDF.at[row,'envfrom'])}
    
    
    
    CountHMMCompDict={x:goodBlocks2[x]['hmmMatchTotal'] for x in goodBlocks2.keys() if len(goodBlocks2[x]['rows'])==2 and (goodBlocks2[x]['hmmMatchTotal'] >15000)}
    newRowHoldings=[]
    if len(CountHMMCompDict)>0:
        
        #newRowHoldings=[]
        for grouping in CountHMMCompDict.keys():
            
            reworkDF = redoDF.loc[goodBlocks2[grouping]['rows']].reset_index().copy()
    
            if list(reworkDF['strand'])[0] == '+' and len(set(reworkDF['strand'])) ==1 and int(reworkDF.iloc[0]['hmm_to']) > int(reworkDF.iloc[1]['hmmfrom']) and int(reworkDF.iloc[1]['hmm_to']) > int(reworkDF.iloc[0]['hmm_to']):
                orientation = list(reworkDF['strand'])[0]
                newRowHoldings.append(['DROPME', str(reworkDF.iloc[0]['target_name']), 'TSPY', int(reworkDF.iloc[0]['hmmfrom']),
                                  int(reworkDF.iloc[1]['hmm_to']), min([int(x) for x in reworkDF['envfrom']]),
                                  max([int(x) for x in reworkDF['env_to']]), orientation, sum([float(x) for x in reworkDF['E_value']]),
                                 sum([float(x) for x in reworkDF['score']]), 'SpecialMerge', '*', '*'])
            
            elif list(reworkDF['strand'])[0] == '-' and len(set(reworkDF['strand'])) ==1 and int(reworkDF.iloc[1]['hmm_to']) > int(reworkDF.iloc[0]['hmmfrom']) and int(reworkDF.iloc[0]['hmm_to']) > int(reworkDF.iloc[1]['hmm_to']):
                orientation = list(reworkDF['strand'])[0]
                newRowHoldings.append(['DROPME', str(reworkDF.iloc[0]['target_name']), 'TSPY', int(reworkDF.iloc[1]['hmmfrom']),
                                  int(reworkDF.iloc[0]['hmm_to']), min([int(x) for x in reworkDF['envfrom']]),
                                  max([int(x) for x in reworkDF['env_to']]), orientation, sum([float(x) for x in reworkDF['E_value']]),
                                 sum([float(x) for x in reworkDF['score']]), 'SpecialMerge', '*', '*'])

        
            else:
                continue
    else:
        pass
    
    if len(newRowHoldings)>0:
        newSpecialDF = pd.DataFrame(data=newRowHoldings).drop(columns=[0])
        # Concatenating DataFrames
        combined_df = pd.concat([refinedDF, newSpecialDF], ignore_index=True)
        combined_df[5]=combined_df[5].astype(int)
        combined_df[6]=combined_df[6].astype(int)
        combined_df.sort_values(by=[1,5,6],inplace=True)
        finalCombinedDF = combined_df.reset_index().drop(columns=['index']).copy()
        finalCombinedDF.to_csv(str(args.output))
    else:
        refinedDF.to_csv(str(args.output))

if __name__=="__main__":
    main()
