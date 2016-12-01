# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 11:49:52 2016

@author: simon
"""

from pymira import amiramesh

class SpatialGraph(amiramesh):
    
    def __init__(self):
        return amiramesh.__init__(self)
        
    def read(self,*args,**kwargs):
        if not amiramesh.read(self,*args,**kwargs):
            return False
        if self.get_parameter_value("ContentType")!="HxSpatialGraph":
            print 'Wanring: File is not an Amira SpatialGraph!'
            return True
    
    def remove_intermediate_nodes(self):
        pass

#PRO REMOVE_INTERMEDIATE_NODES
#
#;dir = 'C:\Users\simon\Dropbox\Mesentery\'
#;graphFile = dir+'Flow2AmiraPressure.am'
#dir = 'C:\Users\simon\Dropbox\Tumour\'
#graphFile = dir+'Flow2AmiraPressure.am'
#ofile = dir+'Flow2Amira_RIN.am'
#
#PRINT, 'Reading graph...'
#amdata = READ_SPATIALGRAPH_AM(graphFile,PARAMETER=param, VALUE=paramVal, COMP=comp, TRANSFORM=transform)
#
#nodeCoords = (*amData[0].val)
#nnode = N_ELEMENTS(nodeCoords[*,0])
#conn = (*amdata[1].val)
#nseg = amdata[1].nval
#edgePoints = (*amdata[3].val)
#nEdgePoint = amdata[3].nval
#
#restore_conn = 1
#f = dir+'nconn.sav'
#IF restore_conn EQ 0 THEN BEGIN
#   PRINT, 'Calculating number of node connections...'
#   nconn = NUMBER_OF_NODE_CONNECTIONS(amdata,CONN=conn)
#   SAVE, nconn,FILE=f
#ENDIF ELSE BEGIN
#  RESTORE, f
#ENDELSE
#
#sInt = WHERE(nconn EQ 2,nInt) ;Nodes with exactly two edges connected
#IF nInt EQ 0 THEN BEGIN
#  PRINT,'No intermediate nodes to remove!
#  RETURN
#ENDIF
#
#;Find intermediate nodes connected to branches or endpoints
#restore_strtconn = 1
#f = dir+'strtconn.sav'
#IF restore_strtconn EQ 0 THEN BEGIN
#  PRINT, 'Identifying intermediate nodes connected to branches or endpoints...'
#  strtConn = LONARR(nInt) - 1
#  FOR i=0,nInt-1 DO BEGIN
#    curNode = sInt[i]
#    ;Find which segments are connected to the current node
#    sSeg = WHERE(conn[*,0] EQ curNode OR conn[*,1] EQ curNode,cSeg)
#    added = 0
#    FOR j=0,cSeg-1 DO BEGIN
#       ;IF sSeg[j] EQ 7720 THEN stop
#       curConn = conn[sSeg[j],*] ;Current segment nodes
#       ;if curconn[0] EQ 7538 AND curconn[1] EQ 7539 then stop
#       ;How many edges are connected to the nodes at the end of the current segment?
#       connBr = nconn[curConn]
#       ;If one or other isn't an intermediate node (2 edges), add it to the list as a starting point
#       IF N_ELEMENTS(connBr) NE 2 THEN STOP
#       stmp = WHERE(connBr NE 2,ctmp)
#       IF ctmp GT 0 THEN BEGIN
#          if sseg[j] EQ 7720 then stop
#          ;print,TRANSPOSE(conn[sseg[j],*]),TRANSPOSE(connbr)
#          strtConn[i] = sSeg[j] ;Store current segment index
#          added += 1
#          ;stmp=where(strtconn eq 7720,ctmp)
#          ;IF ctmp GT 0 THEN STOP
#       ENDIF
#    ENDFOR
#    ;stmp=where(strtconn eq 7720,ctmp)
#    ;IF ctmp GT 0 THEN STOP
#  ENDFOR
#  
#  sStrt = WHERE(strtConn GE 0,nstrt,COMP=sStrtComp)
#  strtConn = strtConn[sStrt]
#  SAVE, strtconn, FILE=f
#ENDIF ELSE BEGIN
#  RESTORE, f
#ENDELSE
#
#converted = LONARR(nInt)
#segVisited = LONARR(nSeg)
#nstrt = N_ELEMENTS(strtconn)
#
#PRINT, 'Number of intermediate points: ',nInt
#PRINT, 'Number of starting points: ',nStrt
#
#keepNode = LONARR(nnode) + 1
#keepEdgePoint = LONARR(nEdgePoint) + 1
#keepSeg = LONARR(nseg) + 1
#
#PRINT, 'Loading a copy of the graph...'
#amDataCopy = READ_SPATIALGRAPH_AM(graphFile)
#
#nodeCoordsCopy = (*amDataCopy[0].val)
#connCopy = (*amDataCopy[1].val)
#nedgePointsCopy = (*amDataCopy[2].val)
#edgePointsCopy = (*amDataCopy[3].val)
#radiusCopy = (*amDataCopy[4].val)
#flowCopy = (*amDataCopy[5].val)
#pressureCopy = (*amDataCopy[6].val)
#
#prog = 0
#progStep = 10
#FOR i=0,nstrt-1 DO BEGIN
#  curProg = i*100./FLOAT(nstrt)
#  IF curProg GE prog THEN BEGIN
#    PRINT,'Progress (%): ',prog
#    prog += progStep
#  ENDIF
#  
#  endloop = 0
#  
#  curSeg = strtConn[i]
#  curSegConn = PARE(conn[curSeg,*])
#  curSegNConn = nconn[curSegConn]
#  strtSegNConn = curSegNConn
#  
#  sStrt = WHERE(curSegNConn NE 2,cStrt,COMP=sNext)
#  IF cStrt NE 1 THEN BEGIN
#    endloop = 1
#    PRINT, 'Error: Started loop in non-terminating intermediate node...'
#    stop
#  ENDIF
#  curEdge = GET_EDGE(amdata,curSeg)
#  edge = [curEdge]
#  segVisited[curSeg] = 1  
#  startNode = PARE(curSegConn[sStrt])
#  curNode = PARE(curSegConn[sNext])
#  
#  sNext = WHERE(curSegConn NE startNode,cNext,COMP=sCur)
#  IF cNext NE 1 THEN STOP
#  revFlag = PARE(sCur[0])
#  
#  keepSeg[curSeg] = 0
#  keepEdgePoint[edge.strtInd:edge.endInd] = 0
#  keepNode[sNext] = 0
#  
#  valid = 0
#  count = 0
#  ;revFlag = []
#  REPEAT BEGIN
#     sLoc = WHERE(conn[*,0] EQ curNode OR conn[*,1] EQ curNode AND segvisited EQ 0,cLoc)
#     IF cLoc EQ 0 THEN BEGIN
#        endloop = 1
#        BREAK
#     ENDIF
#     IF endloop EQ 0 THEN BEGIN
#       IF cLoc GT 1 THEN STOP
#       curSeg = sLoc[0]
#       curSegConn = conn[curSeg,*]
#       curSegNConn = nconn[curSegConn]
#       
#       sNext = WHERE(curSegConn NE curNode,cNext,COMP=sCur)
#       IF cNext NE 1 THEN STOP
#       revFlag = [revFlag,PARE(sCur[0])]
#       nextNode = PARE(curSegConn[sNext])
#       
#       curEdge = GET_EDGE(amdata,curSeg)
#       edge = [edge,curEdge]
#       keepSeg[curSeg] = 0
#       keepEdgePoint[curEdge.strtInd:curEdge.endInd] = 0
#       
#       CASE 1 OF
#          curSegNConn[0] EQ 2 AND curSegNConn[1] EQ 2: BEGIN
#            segVisited[curSeg] = 1
#            keepNode[sCur] = 0
#            curNode = PARE(nextNode)
#          END
#          ELSE: BEGIN
#            ;IF i EQ 70 THEN STOP
#            ;IF revFlag[-1] EQ 1 THEN BEGIN
#            ;  endNode = PARE(curNode)
#            ;ENDIF ELSE BEGIN
#              endNode = PARE(nextNode)
#            ;ENDELSE
#            endloop = 1
#            valid = 1
#          END
#       ENDCASE
#       count += 1
#     ENDIF
#  ENDREP UNTIL endloop EQ 1
#  
#  IF valid EQ 1 THEN BEGIN
#    ;Add new edge into graph
#    IF revFlag[0] EQ 0 THEN BEGIN
#      radius = (*edge[0].radius)
#      flag = (*edge[0].flag)
#      flow = (*edge[0].flow)
#      pressure = (*edge[0].pressure)
#      points = (*edge[0].points)
#    ENDIF ELSE BEGIN
#      radius = REVERSE(*edge[0].radius)
#      flag = REVERSE(*edge[0].flag)
#      flow = REVERSE(*edge[0].flow)
#      pressure = REVERSE(*edge[0].pressure)
#      points = REVERSE(*edge[0].points)
#    ENDELSE
#    
#    FOR j=1,N_ELEMENTS(edge)-1 DO BEGIN
#      IF revFlag[j] EQ 0 THEN BEGIN
#         points = [points,(*edge[j].points)[1:*,*]]
#         radius = [radius,(*edge[j].radius)[1:*]]
#         flag = [flag,(*edge[j].flag)[1:*]]
#         flow = [flow,(*edge[j].flow)[1:*]]
#         pressure = [pressure,(*edge[j].pressure)[1:*]]
#      ENDIF ELSE BEGIN
#         ncur = N_ELEMENTS((*edge[j].radius))
#         points = [points,REVERSE((*edge[j].points)[0:ncur-2,*])]
#         radius = [radius,REVERSE((*edge[j].radius)[0:ncur-2])]
#         flag = [flag,REVERSE((*edge[j].flag)[0:ncur-2])]
#         flow = [flow,REVERSE((*edge[j].flow)[0:ncur-2])]
#         pressure = [pressure,REVERSE((*edge[j].pressure)[0:ncur-2])]
#      ENDELSE
#    ENDFOR
#    
#    stmp1 = WHERE(PARE(points[0,*]) NE nodeCoords[startNode,*],ctmp1)
#    IF ctmp1 GT 0 THEN STOP
#    nel = N_ELEMENTS(points[*,0])
#    stmp2 = WHERE(PARE(points[nel-1,*]) NE nodeCoords[endNode,*],ctmp2)
#    IF ctmp2 GT 0 THEN STOP
#
#    connCopy = [[connCopy],TRANSPOSE([startNode,endNode])]
#    edgePointsCopy = [edgePointsCopy,points]
#    nedgePointsCopy = [nedgePointsCopy,N_ELEMENTS(radius)]
#    radiusCopy = [radiusCopy,radius]
#    flowCopy = [flowCopy,flow]
#    pressureCopy = [pressureCopy,pressure]
#  ENDIF
#ENDFOR
#
#nSegAdded = N_ELEMENTS(conncopy[*,0]) - nSeg
#
#nedgepoints = N_ELEMENTS(radiusCopy)
#amdataCopy[0].val = PTR_NEW(nodeCoordsCopy)
#amdataCopy[0].nval = N_ELEMENTS(nodeCoordsCopy[*,0])
#amdataCopy[1].val = PTR_NEW(connCopy)
#amdataCopy[1].nval = N_ELEMENTS(connCopy[*,0])
#amdataCopy[2].val = PTR_NEW(nedgePointsCopy)
#amdataCopy[2].nval = N_ELEMENTS(nedgePointsCopy)
#amdataCopy[3].val = PTR_NEW(edgePointsCopy)
#amdataCopy[3].nval = nedgepoints
#amdataCopy[4].val = PTR_NEW(radiusCopy)
#amdataCopy[4].nval = nedgepoints
#amdataCopy[5].val = PTR_NEW(flowCopy)
#amdataCopy[5].nval = nedgepoints
#amdataCopy[6].val = PTR_NEW(pressureCopy)
#amdataCopy[6].nval = nedgepoints
#
#SAVE, /VAR, FILE=dir+'RIN_vars.sav'
#
#keepSegNew = [keepSeg,INTARR(nSegAdded)+1]
#
#;skeep = WHERE(keepSegNew EQ 1,cKeep,COMP=sDel,NCOMP=nDel)
#;nSegTotal = nsegAdded+nSeg
#;newSegMapping = INTARR(2,nSegTotal)
#;newSegMapping[0,*] = INDGEN(nSegTotal)
#;newSegMapping[1,sDel] = -1
#;newSegMapping[1,sKeep] = INDGEN(cKeep)
#;
#;newNodeMapping = INTARR(2,nnode)
#;newNodeMapping[0,*] = INDGEN(nnode)
#;newNodeMapping[1,*] = INDGEN(nnode)
#;
#;segCount = 0
#;FOR i=0,nSegTotal-1 DO BEGIN
#;  segNodes_oldInds = PARE(conncopy[skeep[i],*])
#;  IF keepSegNew[i] EQ 0 THEN BEGIN
#;    e = GET_EDGE(amdataCopy,i)
#;    
#;  ENDIF
#;ENDFOR
#;stop
#
#DELETE_EDGE, amdataCopy, keepSegNew
#REMOVE_ORPHAN_NODES, amdatacopy
#
#WRITE_SPATIALGRAPH_AM_2, amdataCopy, FILE=ofile, PARAM=param, VALUE=paramVal, TRANSFORM=transform
#
#END
