# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 14:40:52 2016

@author: simon
"""

import numpy as np
import re

class AmiraMesh(object):
    
    def __init__(self):
        self.parameters = None
        self.paramText = None
        self.fileType = None
        self.definitions = None
        self.fields = None
        self.fieldOffset = None
        self.data = []
        
    def _dataline(self,curLine):
        endPat = '1234567890@ \n'
        chk = len([x for x in curLine if x in endPat])==len(curLine) and '@' in curLine
        marker = None
        markerIndex = 0
        if chk:
            mtch = re.search(r'@\w',curLine)
            if mtch is not None:
                try:
                    marker = mtch.group(0)
                    markerIndex = int(marker.replace('@',''))
                except Exception,e:
                    print repr(e)
                    return chk,marker,markerIndex
        return chk,marker,markerIndex
        
    def _read_file_data(self,content,strt,fin,curDataMarker):
        
        # Get corresponding field and definition entries
        curField = [x for x in self.fields if x['marker']==curDataMarker][0]
        curDef = [x for x in self.definitions if x['name']==curField['definition']][0]
        
        # Update field entry to include data size and shape
        curField['nentries'] = [long(x) for x in curDef['size']]
        curShape = [long(x) for x in curDef['size']]
        if curField['nelements']>1:
            curShape.append(curField['nelements'])
        curField['shape'] = curShape
        
        # Data type
        if 'ascii' in self.fileType.lower():
            if curField['type']=='float':
                dtype = np.dtype('f')
            elif curField['type']=='double':
                dtype = np.dtype('d')
            elif curField['type']=='byte':
                dtype = np.dtype('b')                     
            elif curField['type']=='int':
                dtype = np.dtype('i')
            elif curField['type']=='long':
                dtype = np.dtype('i')
            else: # Default to float
                dtype = np.dtype('f')
        else:
            print('Non-ASCII files not currently supported!')
            return False
            
        # Grab data

        # Check size
        if (fin-strt)!=int(curDef['size'][0]):
            import pdb
            pdb.set_trace()
            print 'Error: Data size does not match header definition! Marker: '+curDataMarker
            return None
        curData = [x.rstrip().split(' ') for x in content[strt:fin] if x]
            
        # Convert to numpy array
        curData = np.array(curData,dtype=dtype)
        # Reshape data
        curData = np.reshape(curData,curField['shape'])
        # Add to field entry
        curField['data'] = curData
        
        return True
        
    def read(self,filename):
        
        with open(filename) as f:
            content = f.readlines()
            
        self.parameters = []
        self.definitions = []
        self.fileType = None
        self.paramText = []
        self.fields = []
        self.fieldOffset = []
        self.data = []
            
        # Remove line returns
        content = [x.strip('\n') for x in content]

        inHeader = True
        inParam = False
        curLine = ''
        firstLineFound = False
        nonEmptyLineCount = 0
        dataFieldCount = 0
        i = -1
        
        for i,curLine in enumerate(content):
            #i += 1 # Use this instead of enumerate
            
            if curLine=='' or curLine==' ':
                emptyLine = True
            else:
                nonEmptyLineCount += 1
                emptyLine = False

            if inHeader:
                # Check for line definition symbols (if line is short)
                if len(curLine)<200:
                    defChk = 'define' in curLine.lower()
                    hashChk = '#' in curLine
                    atChk = '@' in curLine
                    paramChk = 'parameters' in curLine.lower()
                    dataSectionChk,curDataMarker,curMarkerIndex = self._dataline(curLine)
                else:
                    defChk = False
                    hashChk = False
                    atChk = False
                    paramChk = False
                    dataSectionChk = False
                    
                if nonEmptyLineCount==1 and firstLineFound is False:
                    # Look for magic
                    firstLineFound = True
                    ptn = re.compile("\s*#\s*amiramesh")
                    magicMtch = ptn.match(curLine.lower()) is not None
                    if not magicMtch:
                        print('Error: Not a supported Amira file!')
                        return False
                    spl = curLine.split('AmiraMesh')
                    self.fileType = ''.join(spl[1:]).strip()
                    
                # Parameter line
                if paramChk:
                    self.parameters = []
                    inParam = True
                    bracketCount = 0
                    
                if inParam:
                    self.paramText.append(curLine)
                    spl = str.split(curLine,' ')
                    spl = [x for x in spl if x!='']
                    nspl = len(spl)
                    if '{' in curLine:
                        bracketCount += 1
                    if '}' in curLine:
                        bracketCount -= 1
                    #print 'Bracket: ',bracketCount
                    if nspl>1 and '{' not in curLine and '}' not in curLine: # Nested parameter structure not yet supported
                    
                        curParam = spl[0]
                        if '"' in curLine:
                            val = ' '.join(spl[1:])
                            val = val.replace('"','').rstrip()
                            self.parameters.append({'parameter':curParam,'value':val})
                        elif curParam=='ContentType':
                            spl2 = str.split(spl[1],'\"')
                            curParamVal = spl2[0]
                            self.parameters.append({'parameter':curParam,'value':curParamVal})                        
                        elif curParam=='TransformationMatrix':
                            matr = [float(re.sub(r'\W+','',x)) for x in spl[1:]]
                            self.parameters.append({'parameter':curParam,'value':matr})                        
                        elif curParam=='BoundingBox':
                            bbox = [float(re.sub(r'\W+','',x)) for x in spl[1:]]
                            self.parameters.append({'parameter':curParam,'value':bbox})                        
                        else:
                            val = spl[1:]
                            self.parameters.append({'parameter':curParam,'value':val})
    
                    if bracketCount==0:
                        inParam = False
                            
                if defChk and not inParam: # Definition line
                    spl = str.split(curLine, ' ')
                    nspl = len(spl)
                    if nspl>2:
                        self.definitions.append({'name':spl[1],'size':[long(x) for x in spl[2:]]})
                    else:
                        print('Error: Unexpected formatting on definition line!')
                        return False
                        
                if atChk and not inParam and not dataSectionChk: # Field line
                    spl = curLine.split(' ')
                    spl = [x for x in spl if x!='']
                    
                    defMtch = [x for x in self.definitions if x['name']==spl[0]]
                    
                    if len(defMtch)>0:
                        spl = re.split(r"[{}]+", curLine)
                        spl = [x for x in spl if x!='']
                        fieldDef = spl[0].strip()
                        fieldType = spl[1].strip()
                        typeSpl = spl[1].split(' ')
                        typeSpl = [x for x in typeSpl if x!='']
                        if len(typeSpl)==2:
                            fieldName = typeSpl[1].strip()
                            fieldDataType = typeSpl[0].strip()
                            m = re.split(r"\[(\w+)\]", fieldDataType)
                            m = [x for x in m if x!='']
                            if len(m)==2:
                                fieldDataType = m[0].strip()
                                nel = int(m[1].strip())
                            else:
                                nel = 1
                        else:
                            print('Error: Unexpected formatting on variable definition line! ',curLine)
                            return None
                        spl2 = re.split(r"[()]+", spl[2])
                        spl2 = [x for x in spl2 if x!='']
                        marker = spl2[0].strip()
                        if len(spl2)>1:
                            fieldTypeInfo = spl2[1].split(',')
                            fieldTypeInfo = [x for x in fieldTypeInfo if x!='']
                            if len(fieldTypeInfo)!=2:
                                import pdb
                                pdb.set_trace()
                        else:
                            fieldTypeInfo = [None,None]
                         
                        self.fields.append({'name':fieldName,'marker':marker, \
                                             'definition':fieldDef,'type':fieldDataType, \
                                             'encoding':fieldTypeInfo[0],'length':fieldTypeInfo[1],\
                                             'nelements':nel,'nentries':None})
                    else:
                        pass
                
                if dataSectionChk: # Come out of header and move on to data
                    self.fieldOffset.append(i)         
                    inHeader = False
                    
            else:
                # The first time a data section marker is found (usually '@1') is during
                # reading the header (it's hopw we know the header has ended). When the next
                # data section marker (usually '@2') is found, we should end up here.
                # This means that we are always sorting through the data once we've been
                # through it already.
                dataSectionChk,newDataMarker,newMarkerIndex = self._dataline(curLine)
                eofChk = i==(len(content)-1) # Check for end of file
                
                if dataSectionChk or eofChk:
                    dataFieldCount += 1
                    self.fieldOffset.append(i)
                    strt = self.fieldOffset[dataFieldCount-1]+1
                    fin = self.fieldOffset[dataFieldCount]-1
                    if eofChk:
                        fin += 1
                    print curDataMarker
                    if not self._read_file_data(content,strt,fin,curDataMarker):
                        return False
                    
                    curDataMarker = newDataMarker # Increment the marker
        
        return True
        
    def test_read(self):
        am_file = r'G:\OPT\19.09.2016 Subcuts Tom\LSM2\lsm2-files\LSM2_rec2.labels'
        # Spatial graph
        am_file = r'G:\OPT\19.09.2016 Subcuts Tom\LSM2\LSM2_bgrem_frangi_response_skel.SptGraph.am'
        tmp = self.read(am_file)
        import pdb
        pdb.set_trace()