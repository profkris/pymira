path = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\\LS#4\\'
f = path+'ls4_vessel_seg.Smt.SptGraph.am'
ps = 8.21
ofile = path + 'LS4_spatialGraph_scaled.am'

path = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\\LS#1\\'
f = path+'LS1_spatialGraph_scaled.am'
ps = None
ofile = None

# path = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\\SW#1\\'
# f = path+'SW1_vessel_seg.Smt.SptGraph.am'
# ps = 4.87
# ofile = 'SW1_spatialGraph_scaled.am'

# path = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\\SW#2\\'
# f = path+'SW2_vessel_seg.Smt.SptGraph.am'
# ps = 6.45
# ofile = 'SW2_spatialGraph_scaled.am'

# path = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\\SW#3\\'
# f = path+'SW3_vessel_seg_skel.SptGraph_with_radius.am'
# ps = 11.02
# ofile = 'SW3_spatialGraph_scaled.am'

from pymira import spatialgraph, statistics
graph = spatialgraph.SpatialGraph()
print('Reading graph... {}'.format(f))
graph.read(f)
print('Graph read')
if ps is not None:
    graph.rescale_coordinates(ps,ps,ps)
    graph.rescale_radius(ps)
    graph.write(ofile)
stats = statistics.Statistics(graph)
stats.do_stats(output_directory=path)