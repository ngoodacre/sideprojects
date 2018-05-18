wdir='C:\\Users\\Norman\\My Documents\\VCU_research\\yeast_eDUFs'
## vector containing gene names, corresponding to G 
#G_names=[,,]

## array containing domain names, corresponding to D
#D_names=[[],[],[]]

## array containing domain names, corresponding to Dm
#Dm_names=[,,,]

## array containing essentiality information for each domain
#Dm=[,,]

## S provides, in inverse order, the indices of all instances of domain k,
## ... for all m domains in Dm
#S=[]


#G is a vector containing essentiality information for all genes (0=non-essential,1=essential)
#D is an array containing essentiality information for all domains, where rows (index level 1)
#.. correspond to genes, and within a row, indices correspond to domains (0=non-essential,1=essential)
def estimate_D(Gcf,Gcf_corr,D,Dm,i,pfam,sigma_k,FER,FNR):
    #print 'estimating D'
    g_i1,g_i2=Gcf[i]
    d_j=Dm[pfam]
    if d_j==0:
        return "non-essential domain"
    #k=Dm_nameindex[Dname]
    #D_given_param=likelihood_D_given_param(Dm,k,sigma)
    G_given_D_andparam=likelihood_G_given_D_andparam(Gcf,Gcf_corr,i,D,Dm,FER,FNR)
    G_given_param=likelihood_G_given_param(Gcf,i,FER,FNR)
    p_ij=sigma_k*G_given_D_andparam/G_given_param
    #p_ij=sigma_k*G_given_D_andparam
    return p_ij

def maximize_param(Gcf,Gcf_corr,D,Dm,Dm_names,S,sigma,FER,FNR):
    print 'maximizing param'
    new_sigma=[]
    print len(sigma)
    print len(S)
    for k,sigma_k in enumerate(sigma):
        print k
        Sk=S[k]
        size_k=len(Sk)
        p_k=float(0)
        pfam=Dm_names[k]
        for i in Sk:
            p_di=estimate_D(Gcf,Gcf_corr,D,Dm,i,pfam,sigma_k,FER,FNR)
            if p_di=="non-essential domain":
                continue
##            if p_di>1:
##                print pfam
##                print p_di
            p_k+=p_di
        try:
            new_sigma_k=p_k/size_k
        except ZeroDivisionError:
            new_sigma_k=float(0)
        new_sigma.append(new_sigma_k)
    return new_sigma

# updates the domain essentiality vector (Dm), using new sigma parameter from most recent maximization step
def update_Dm(Dm,Dm_names,new_sigma,sigma_threshold):
    new_Dm=dict()
    for k,Dm_name in enumerate(Dm_names):
        sig=new_sigma[k]
        if sig>sigma_threshold:
            new_Dm[Dm_name]=1
        else:
            new_Dm[Dm_name]=0
    return [Dm,new_Dm]

# Dm is a vector containing essentiality information for all domains (0=non-essential,1=essential)
# output is L(D|param)
def likelihood_D_given_param(Dm,k,sigma):
    print 'calculating P(D|param)'
    dm_k = Dm[k]
    if dm_k==1:
        L_k=sigma
    if dm_k==0:
        L_k=1-sigma
    return L_k

# "i" is the gene index
# output is L(G|D,param) for a particular g_i
def likelihood_G_given_D_andparam(Gcf,Gcf_corr,i,D,Dm,FER,FNR):
    #print 'calculating P(G|D,param)'
    g_i1,g_i2=Gcf[i]
    g_essent=Gcf_corr[i]
    Dnames1=D[g_i1]
    Dnames2=D[g_i2]
    has_essential_d1=False
    has_essential_d2=False
    for dname1 in Dnames1:
        d_essent=Dm[dname1]
        if d_essent==1:
            has_essential_d1=True
    for dname2 in Dnames2:
        if d_essent==1:
            has_essential_d2=True
    if g_essent==0:
        if has_essential_d1 and has_essential_d2:
            L_i=FNR
        else:
            L_i=1-FER
    else:
        if has_essential_d1 and has_essential_d2:
            L_i=1-FNR
        else:
            L_i=FER
    return L_i

def likelihood_G_given_param(Gcf_corr,i,FER,FNR):
    #print 'calculating P(G|param)'
    g_essent=Gcf_corr[i]
    if g_essent==0:
        L_i=1-FNR
    else:
        L_i=1-FER
    return L_i

def calc_FER(Gcf,Gcf_corr,D,Dm):
    print 'calculating FER'
    num=float(0)
    denom=float(0)
    new_FER=float(0)
    for i,g_i in enumerate(Gcf):
        g_i1=g_i[0]
        g_i2=g_i[1]
        g_essent=Gcf_corr[i]
        Dnames1=D[g_i1]
        Dnames2=D[g_i2]
        for ii,Dnames in enumerate([Dnames1,Dnames2]):
            diff_sum_D_ii=float(0)
            g_i=[g_i1,g_i2][ii]
            for j,dname in enumerate(Dnames):
                d_j=Dm[dname]
                num+=(1-g_essent)*(1-d_j)
                denom+=(1-d_j)
##            num+=(1-g_essent)*(1-diff_sum_D_ii)
##            denom+=(1-diff_sum_D_ii)
    try:
        new_FER=num/denom
    except ZeroDivisionError:
        new_FER=float(0)
    return new_FER 

def calc_FNR(Gcf,Gcf_corr,D,Dm):
    print 'calculating FNR'
    num=float(0)
    denom=float(0)
    for i,g_i in enumerate(Gcf):
        g_i1=g_i[0]
        g_i2=g_i[1]
        g_essent=Gcf_corr[i]
        Dnames1=D[g_i1]
        Dnames2=D[g_i2]
        for ii,Dnames in enumerate([Dnames1,Dnames2]):
            diff_sum_D_ii=float(0)
            g_i=[g_i1,g_i2][ii]
            for j,dname in enumerate(Dnames):
                d_j=Dm[dname]
                num+=g_essent*(1-d_j)
                denom+=(1-d_j)
##            num+=g_essent*diff_sum_D_ii
##            denom+=diff_sum_D_ii
    try:
        new_FNR=num/denom
    except ZeroDivisionError:
        new_FNR=float(0)
    return new_FNR
    
