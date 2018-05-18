

deg_dir='C:\\Users\\Norman\\Documents\\VCU_research\\DEG'
deg_annotation_filenames=[deg_dir+'\\'+'degannotation-'+physymbol+'.dat' for physymbol in ['e','a','p']]


def parse_DEG_annotation_file(deg_annotation_filename,infotypes):
    inf=open(deg_annotation_filename)
    gis=[]
    for i,line in enumerate(inf):
        if i==0:
            continue
        sl=line.strip().split('\t')
        acc=sl[2]
        ##### extract basic info
        if acc.startswith('GI'):
            acc=acc.split(':')[1]
        gis.append(acc)
    return gis

infotypes=[]
gis_e=parse_DEG_annotation_file(deg_annotation_filenames[0],infotypes)
gis_a=parse_DEG_annotation_file(deg_annotation_filenames[1],infotypes)
gis_p=parse_DEG_annotation_file(deg_annotation_filenames[2],infotypes)
