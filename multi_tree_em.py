from cogent import LoadTree, LoadSeqs, DNA
from cogent.evolve.models import JC69
from simulator import simulate_alignment
#import matplotlib.pyplot as plt
import numpy as np
import copy

##input is a list that contains N elements each is the loglikelihood value that the site belong to the group
##for example, when N=3, there will be three values, it is one site belong to branch-length set 1,2,3's loglikelihood value

###output the group that this site belongs to

###process:convert to likelihood, sum up and normalise (likelihood/sum) to make all the value add up to 1
###generate a random number, if normalised value is 0.1,0.8,0.1,then if random number<0.1, return 1, 0.1-0.9 return 2, otherwise return 3

def normalize_and_decide(prob):
    exp_probs=[]
    normaliza_prob=[0]*len(prob);
    sum_of_previous_prob=[0]*len(prob);
    for i in range(len(prob)):
        exp_prob=np.exp(prob[i])
        exp_probs.append(exp_prob)
    _sum = sum(exp_probs)
    for i in range(len(exp_probs)):
        normaliza_prob[i]=exp_probs[i]/_sum
        sum_of_previous_prob[i]=sum(normaliza_prob);
    ran_num=np.random.random();
    for i in range(len(prob)):
        if sum_of_previous_prob[i]>ran_num:
            return i+1;


    return 0;

####input the 2000 sites and N trees, decide which site belong to which tree

def expectation_singlesite(aln, trees):
    modle = JC69()
    result = []
    ####contain len(aln) element, each is a list single_site_ll
    #####contain n element, from 0---n-1 each indicates the loglikelihood to each tree

    n=len(trees)
    for i in range(len(aln)):
        single_site_ll = []
        for j in range(n):
            lf=modle.makeLikelihoodFunction(trees[j]);
            lf.setAlignment(aln[i]);
            prob=lf.getLogLikelihood()
            single_site_ll.append(prob)
        single_result=normalize_and_decide(single_site_ll)
        result.append(single_result)
    print result
    return result

###input site assignment result and aln, the branch lengths N set
###return the optimised likelihood and optimised n branch lengths set value

def optimization(result, aln, trees):

    # get the sites for each tree according to the assignments
    known_tree = '((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
    n=len(trees);
    aln_single=[]
    for i in range(n):
        aln1 = LoadSeqs(data=[('a', ''), ('c', ''),('b', ''), ('d', '')], moltype=DNA)
        aln_single.append(aln1);
    for i in range(len(aln)):
        tree_index = result[i];
        if(tree_index!=0):
            aln_single[tree_index-1]=aln_single[tree_index-1]+aln[i]
    modle = JC69()
    #print aln_single
    # calculate the likelihood and do optimization. optimise will generates
    # new tree parameters
    likelihood_single_tree=[]
    trees_param=[]
    for j in range(n):
        single_tree_param=[]
        lf1 = modle.makeLikelihoodFunction(trees[j])
        lf1.setAlignment(aln_single[j])
        lf1.optimise(local=True)
            #print lf1
        likelihood1 = lf1.getLogLikelihood()
        likelihood_single_tree.append(likelihood1)
         
        single_tree_param.append((lf1.getParamValue('length', 'a') + lf1.getParamValue('length', 'c')) / 2.0)
        single_tree_param.append((lf1.getParamValue('length', 'b') + lf1.getParamValue('length', 'd')) / 2.0)
        single_tree_param.append(lf1.getParamValue('length', 'edge.1') + lf1.getParamValue('length', 'edge.0'))
            #print tree_param
            #trees[j]= LoadTree(treestring=known_tree % (tree_param[0],tree_param[1],tree_param[2] / 2.0,tree_param[0],tree_param[1],tree_param[2]/2.0))
            #print trees[j]
        trees_param.append(single_tree_param)

    likelihood=sum(likelihood_single_tree)
    #print likelihood_single_tree
    
    return trees_param, likelihood

###input aln and number of tree sets
###output is the best optimised parameters, best site assignment result and largest likelihood for N set

def multi_tree(n,aln):
    best_ll = -np.inf
    best_result=[]
    known_tree = '((a:%f, b:%f):%f,(c:%f,d:%f):%f);'
    tree = LoadTree(treestring='((a,b),(c,d))')
    tree_param=[];
    ###generate 3*n random numbers as initial result
    for i in range(n):
        single_tree_param=[]
        for m in range(3):
            single_tree_param.append(np.random.exponential(-1.0 / np.log(0.05)));
        tree_param.append(single_tree_param);

   ###iteration between E and M
    for i in range(200):
        print i
        trees=[]
        for j in range(n):
            single_tree= LoadTree(treestring=known_tree % (tree_param[j][0],tree_param[j][1],tree_param[j][2] / 2.0,tree_param[j][0],tree_param[j][1],tree_param[j][2]/2.0))
            trees.append(single_tree)
        result = expectation_singlesite(aln, trees)
        optimization_result = optimization(result, aln, trees)
        tree_param = optimization_result[0]
        ###if likelihood is larger than the current largest one, keep its site assignment result and branch lengths result
        if optimization_result[1] > best_ll:
            best_ll = optimization_result[1]
            best_param = copy.copy(tree_param)
            best_result=result
            print best_result
            print best_ll
            print best_param

    return best_ll,best_param,best_result

def one_experiment():
    ############simulate sites that under three sets of branch lengths
    ###simulate a 2000bp alignment p1=0.19 q1=0.05,r1=r2=0.1,p2=0.01 q2=0.95
    ############################## p3=0.19 q3=0.02 r3=0.8 p4=0.01 q4=0,38 r4=0.8
    aln1, tree1 = simulate_alignment(0.1, 0.5, 0.1, 1000, 0.9, 'fixed')
    aln2, tree1 = simulate_alignment(0.1, 0.2, 0.8, 1000, 0.9, 'fixed')
    aln=aln1+aln2
    likelihoods=[]
    bic=[]
    branch_lengths=[]
    results=[]
    ### N from 2-6
    for i in range(2,7):
        ll,tree_param,result=multi_tree(i,aln);
        print ll,tree_param,result
        bic_n=np.log(2000)*i*3-2*ll
        likelihoods.append(ll)
        bic.append(bic_n)
        branch_lengths.append(tree_param)
        results.append(result)
    return bic,likelihoods,branch_lengths,results


bic,likelihoods,branch_lengths,results=one_experiment();
print bic
print likelihoods
print branch_lengths
print result
np.savetxt('bic.txt', bic)
np.savetxt('likelihoods.txt', likelihoods)
np.savetxt('result.txt', results)
np.savetxt('branch_lenghths.txt', branch_lengths)





























