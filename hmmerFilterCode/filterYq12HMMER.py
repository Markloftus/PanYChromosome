import pandas as pd
import numpy as np
import os
import collections
from collections import defaultdict
import gzip
import argparse as ap

def merge_pairs(pairs):
    merged_pairs = []
    current_pair = pairs[0]

    for i in range(1, len(pairs)):
        if pairs[i][2] == current_pair[2] and pairs[i][3] == current_pair[3] and pairs[i][0] == current_pair[1]:
            current_pair[1] = pairs[i][1]
        else:
            merged_pairs.append(current_pair)
            current_pair = pairs[i]

    merged_pairs.append(current_pair)
    return merged_pairs



def main():
    global args
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, help='What is the path to the input folder /my/input/folder/hmmerFiles/')
    parser.add_argument("-o", "--output", dest='output', required=True, help='Give the full path to output the dataframe')
    args = parser.parse_args()
    finalBlocks=[]
    # File path to nhmmer output
    tblout_folder = str(args.input)
    
    goodYq12Files=['Yqhet_3k1bp.hmmer-tblout.txt.gz',
    'DYZ18_Yq.hmmer-tblout.txt.gz',
    'DYZ1_Yq.hmmer-tblout.txt.gz',
    'Yqhet_2k7bp.hmmer-tblout.txt.gz',
    'DYZ2_Con.hmmer-tblout.txt.gz'
    ]
    # Dictionary to store the best hits per target coordinate
    
    header=['target_name','accession', 'query_name','accession2', 'hmmfrom', 'hmm_to', 'alifrom','ali_to','envfrom', 'env_to', 'sq_len', 'strand', 'E_value', 'score' ,'bias', 'description_of_target']
    # Parse the nhmmer tblout
    for hmmerFile in os.listdir(tblout_folder):
        if '.hmmer-tblout.txt.gz' in hmmerFile:
            if '.'.join(hmmerFile.split(".")[2:]) in goodYq12Files:
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
    finalDFhalf = finalDF[(finalDF[8]>2) & (finalDF[6]>=0.9)].copy()
    #flag=0
    #flag2=0
    #goodRows=[]
    #for row,element in zip(finalDFhalf.index, finalDFhalf[1]):
    #    if (element) =='Yqhet_2k7bp' and flag==0 and flag2==0:
    #        flag+=1
    #        goodRows.append(row)
    #    elif element =='Yqhet_2k7bp' and flag>0 and flag2>0:
    #        continue
    #    else:
    #        if flag>0 and flag2==0:
    #            flag2+=1
    #            goodRows.append(row)
    #        else:
    #            goodRows.append(row)
    #finalDF2 = finalDFhalf.loc[goodRows]
    goodRows=[]
    for contig in set(finalDFhalf[0]):
        flag=0
        flag2=0
        temporaryDFContig = finalDFhalf[finalDFhalf[0]==contig].copy()
        for row,element in zip(temporaryDFContig.index, temporaryDFContig[1]):
    
            if (element) =='Yqhet_2k7bp' and flag==0 and flag2==0:
                flag+=1
                goodRows.append(row)
            
            elif element =='Yqhet_2k7bp' and flag>0 and flag2>0:
                continue
            else:
                if flag>0 and flag2==0:
                    flag2+=1
                    goodRows.append(row)
                else:
                    goodRows.append(row)
                
    finalDF2 = finalDFhalf.loc[goodRows]

    # Start to do the merging process
    coordinateDict={}
    dropColumns=[]
    coordRowNum=0
    for row in finalDF2.index:
        coordinateDict[row]={'RowNum':coordRowNum, 'Element': str(finalDF2.at[row,1]), 'Contig':str(finalDF2.at[row,0]) ,'Strand':str(finalDF2.at[row,2]), 'Start':int(finalDF2.at[row,3]), 'End':int(finalDF2.at[row,4]),'Drop':'No', 'Trim':'No', 'Rows':'TEMP'}
        coordRowNum+=1
    #Compare coordinates that overlap and do some trimming/dropping
    for row in coordinateDict.keys():
        
        for row2 in coordinateDict.keys():
            
            if row == row2 or coordinateDict[row]['Drop']=='Yes' or coordinateDict[row2]['Drop']=='Yes':
                continue
                
            else:
                
                #If the contig is the same and the elements are different then we need to do something
                if (coordinateDict[row]['Contig'] == coordinateDict[row2]['Contig']) and (coordinateDict[row]['Element'] != coordinateDict[row2]['Element']):
                    
                    #test the range overlap
                    range1 = range(coordinateDict[row]['Start'], coordinateDict[row]['End']+1)
                    range2 = range(coordinateDict[row2]['Start'], coordinateDict[row2]['End']+1)
                    xs = set(range1)
                    overlappingCoordinates = xs.intersection(range2)
    
                    #If the range overlap is over 1000 bases lets do some trimming
                    if len(overlappingCoordinates) > 1000:
    
                        #Figure out which key is more variable
                        minDict={row:float(finalDF2.at[row,7]),row2:float(finalDF2.at[row2,7])}
                        maxKey = max(minDict, key=minDict.get)
                        minKey = min(minDict, key=minDict.get)
    
                        newList=[]
                        if coordinateDict[maxKey]['RowNum'] > coordinateDict[minKey]['RowNum']:
                            for hit in finalDF2.at[maxKey,9]:
                                if int(hit[9]) < int(coordinateDict[minKey]['End']):
                                    continue
                                else:
                                    newList.append(hit)
                        else:
                            for hit in finalDF2.at[maxKey,9]:
                                if int(hit[10]) > int(coordinateDict[minKey]['Start']):
                                    continue
                                else:
                                    newList.append(hit)
    
                        if len(newList)>0:
                            coordinateDict[maxKey]['Trim']='Yes'
                            coordinateDict[maxKey]['Rows']=newList
    
                            newStart = min([int(newHitItem[9]) for newHitItem in newList])
                            newEnd = max([int(newHitItem[9]) for newHitItem in newList])
                            coordinateDict[maxKey]['Start']= newStart
                            coordinateDict[maxKey]['End']= newEnd
    
    
                        else:
                            coordinateDict[maxKey]['Drop']='Yes'
                            dropColumns.append(maxKey)
                    #There is no overlap and just continue to next comparison
                    else:
                        continue
      
                #If the contigs are different or the elements are the same ignore for now
                else:
                    continue
    
    #Build a new dataframe with the remaining coordinates
    finalDF3 = finalDF2.loc[[x for x in finalDF2.index if x not in dropColumns]].copy()
    for row in finalDF3.index:
        if coordinateDict[row]['Trim']=='Yes':
            finalDF3.at[row,3]=coordinateDict[row]['Start']
            finalDF3.at[row,4]=coordinateDict[row]['End']
            finalDF3.at[row,5]= abs(coordinateDict[row]['Start']-coordinateDict[row]['End'])+1
            finalDF3.at[row,8]=len(coordinateDict[row]['Rows'])
            finalDF3.at[row,9]=coordinateDict[row]['Rows']
        else:
            continue
    
    # Make a copy to avoid modifying the original
    result_df = finalDF3.sort_values(by=[0,3]).reset_index().drop(columns=['index']).copy()
    # Initialize list to store merged blocks
    merged_blocks = []
    
    #Find out which adjacent rows can be merged
    for row in result_df.index:
        if row == len(result_df)-1:
            continue
        else:
            rowElementName = result_df.at[row,1]
            contigName = result_df.at[row,0]
            nextRow = row+1
            if rowElementName == result_df.at[nextRow,1] and contigName == result_df.at[nextRow,0]:
                merged_blocks.append([row, nextRow, rowElementName, contigName]) 
            else:
                continue
    mergeTheseRows = merge_pairs(merged_blocks) 
    
    #Merge the data of the rows
    newRows=[]
    usedRows=[]
    for mergedList in mergeTheseRows:
    
        contig = mergedList[3]
        elementName = mergedList[2]
        
        orientations=[]
        rowsInfo =[]
        for number in range(mergedList[0], mergedList[1]+1):
            usedRows.append(number)
            orientations.append(result_df.at[number,2])
    
            for item in result_df.at[number,9]:
                rowsInfo.append(item)
    
        newRowCount= len(rowsInfo)
        
        start = result_df.at[mergedList[0], 3]
        end = result_df.at[mergedList[1], 4]
    
        newTempDF = pd.DataFrame(rowsInfo)
        #make a total length dictionary to see what proportion of the block is zero evalue
        totalLengthDict = {x:0 for x in range(min(newTempDF[9]), max(newTempDF[10])+1)}
        
        for row in newTempDF.index:
            if newTempDF.at[row,13]==0.0:
                for number in range(int(newTempDF.at[row,9]), int(newTempDF.at[row,10])+1):
                    totalLengthDict[number]=1
            else:
                 continue
                        
        totalSize = sum(totalLengthDict.values())/len(totalLengthDict)
        
        #What is the standard deviation of the distances between hits
        blockSTD = np.std(newTempDF.iloc[1:][18])
    
        if len(set(orientations)) ==1:
            newRows.append([contig, elementName, [x for x in orientations][0], start, end, newRowCount, rowsInfo, totalSize, blockSTD])
        else:
            newRows.append([contig, elementName, 'PossibleInversion' , start, end, newRowCount, rowsInfo, totalSize, blockSTD])
    
    #Add back in the rows that didnt have merging issues
    for row in result_df.index:
        if row in usedRows:
            continue
        else:
            newRows.append([result_df.at[row,0],result_df.at[row,1], result_df.at[row,2], result_df.at[row,3],
                            result_df.at[row,4], result_df.at[row,8], result_df.at[row,9], result_df.at[row,6], result_df.at[row,7]])
    
    #Build a final dataframe and resort based on contig and coordinates
    FilteredFinalDF = pd.DataFrame(data=newRows)
    FilteredFinalDF[3] = FilteredFinalDF[3].astype(int)
    FilteredFinalDF[4] = FilteredFinalDF[4].astype(int)
    FilteredFinalDF.sort_values(by=[0,3,4], inplace=True)
    FilteredFinalDF[[0,1,2,3,4,5,7,8]].to_csv(str(args.output))

if __name__=="__main__":
    main()
