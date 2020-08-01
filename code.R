#####均匀设计########
###CD中心化偏差####
CD_function<-function(P){
  n = nrow(P)
  s = ncol(P)
  #线性变换
  P = (P-0.5)/n
  #第一部分
  CD_1 = (13/12)^s
  ###第一部分结束
  #第二部分
  CD_2 = 0
  for(i in 1:n){
    CD_2_multi = 1
    #连乘
    for(j in 1:s){
      #CD_2_multi_iter = 2+abs(P[i,j]-0.5)-abs(P[i,j]-0.5)^2
      CD_2_multi_iter = 1+0.5*abs(P[i,j]-0.5)-0.5*abs(P[i,j]-0.5)^2
      CD_2_multi = CD_2_multi*CD_2_multi_iter
    }
    #print(CD_2_multi)
    CD_2 = CD_2+CD_2_multi
  }
  CD_2 = CD_2*(2/n)
  ###第二部分结束
  #第三部分
  CD_3 = 0
  for(i in 1:n){
    for(k in 1:n){
      CD_3_multi = 1
      #连乘
      for(j in 1:s){
        CD_3_multi_iter = 1+0.5*abs(P[i,j]-0.5)+
          0.5*abs(P[k,j]-0.5)-0.5*abs(P[i,j]-P[k,j])
        CD_3_multi = CD_3_multi*CD_3_multi_iter
      }
      #print(paste(CD_3_multi,'i=',i,'k=',k))
      CD_3 = CD_3+CD_3_multi
    }
    CD_3 = CD_3+CD_3_multi
  }
  CD_3 = CD_3/n^2
  ####第三部分结束###
  #求和开方
  CD = sqrt(CD_1-CD_2+CD_3)
  return(CD)
}

######################
###WD可卷偏差####
WD_function<-function(P){
  n = nrow(P)
  s = ncol(P)
  #线性变换
  P = (P-0.5)/n
  ##第一部分
  WD_1 = (4/3)^s
  ##第一部分结束
  #第二部分
  WD_2 = (1/n)*(3/2)^s
  #第二部分结束
  #第三部分
  WD_3 = 0
  for(i in 1:(n-1)){
    for(k in (i+1):n){
      WD_3_multi = 1
      for(j in 1:s){
        WD_3_multi_iter = (3/2)-abs(P[i,j]-P[k,j])+abs(P[i,j]-P[k,j])^2
        WD_3_multi = WD_3_multi*WD_3_multi_iter                                                   
      }
      WD_3 = WD_3+WD_3_multi
    }
    WD_3 = WD_3+WD_3_multi
  }
  WD_3 = WD_3*2/n^2
  ##第三部分结束
  #求和开放
  WD = sqrt(-WD_1+WD_2+WD_3)
  return(WD)
}

##########
###MD混合偏差####
MD_function<-function(P){
  n = nrow(P)
  s = ncol(P)
  #线性变换
  P = (P-0.5)/n
  #第一部分
  MD_1 = (19/12)^s
  ###第一部分结束
  #第二部分
  MD_2 = 0
  for(i in 1:n){
    MD_2_multi = 1
    #连乘
    for(j in 1:s){
      #MD_2_multi_iter = 2+abs(P[i,j]-0.5)-abs(P[i,j]-0.5)^2
      MD_2_multi_iter = 5/3-0.25*abs(P[i,j]-0.5)-0.25*abs(P[i,j]-0.5)^2
      MD_2_multi = MD_2_multi*MD_2_multi_iter
    }
    #print(MD_2_multi)
    MD_2 = MD_2+MD_2_multi
  }
  MD_2 = MD_2*(2/n)
  ###第二部分结束
  #第三部分
  MD_3 = 0
  for(i in 1:n){
    for(k in 1:n){
      MD_3_multi = 1
      #连乘
      for(j in 1:s){
        MD_3_multi_iter = 15/8-0.25*abs(P[i,j]-0.5)-0.25*abs(P[k,j]-0.5)-
          0.75*abs(P[i,j]-P[k,j])+0.5*abs(P[i,j]-P[k,j])^2
        MD_3_multi = MD_3_multi*MD_3_multi_iter
      }
      #print(paste(MD_3_multi,'i=',i,'k=',k))
      MD_3 = MD_3+MD_3_multi
    }
    MD_3 = MD_3+MD_3_multi
  }
  MD_3 = MD_3/n^2
  ####第三部分结束###
  #求和开方
  MD = sqrt(MD_1-MD_2+MD_3)
  return(MD)
}

##############
####定义同余函数###
#a是被模数，n是模，两者都是正整数，a可以为序列
mod<-function(a,n){
  m = length(a)
  u = rep(0,m)
  for(t in 1:m){
    a_star = a[t]
    if(a_star > n){
      p = floor(a_star/n)
      for(i in 1:p){
        p_star = a_star-i*n
        if(p_star>0 & p_star<n){
          #print(p_star)
          u[t] = p_star
        }
      }
    }else{
      #print(a)
      u[t] = a_star
    }
  }
  #保证余数落入[1，n]的正整数区间
  u[which(u == 0)] = n
  return(u)
}
mod(c(3,10,77,88),10)


##组合数计算公式C^a_b,要求a>b
C<-function(a,b){
  c = a-(b-1)
  #计算分子
  numerator = 1
  for(i in a:c){
    numerator = numerator*i
  }
  #计算分母
  denominator = 1
  for(j in 1:b){
    denominator = denominator*j
  }
  return(numerator/denominator)
}
C(10,2)
#和choose()函数一样，不要了###

#####################################
########好格子点法####
##给定n和s####
n = 15;s = 2
glpm = function(n,s){
  library(schoolmath)
  #初始化均匀设计矩阵
  U = matrix(0,nrow = n,ncol = s)
  U[,1] = c(1:n)
  #构造互质的备选正整数集
  H_n = list()
  i = 1
  for(h in 1:n){
    if(gcd(n,h) == 1){
      H_n[i] = h
      i = i+1
    }
  }
  H_n = unlist(H_n)
  m = length(H_n)
  #如果m >= s则可以继续，否则不能继续
  if(m >= s){
    #number表示一共有多少选法
    #利用choose组合数函数
    number = choose((m-1),(s-1))
    #初始化偏差矩阵,共可产生C^(m-1)_(s-1)个设计矩阵
    CD_mat = matrix(0,nrow = 1,ncol = number)
    WD_mat = matrix(0,nrow = 1,ncol = number)
    MD_mat = matrix(0,nrow = 1,ncol = number)
    #构造生成向量矩阵
    U_gene_mat = matrix(0,nrow = number,ncol = s)
    U_gene_mat[,1] = 1
    #使用combn排列组合函数
    U_gene_mat[,2:s] = t(combn(setdiff(H_n,1),s-1))
    #构造设计矩阵
    for(i in 1:number){
      U[1,]= U_gene_mat[i,]
      for(j in 2:s){
        U[2:n,j] = mod(U[2:n,1]*U[1,j],n)
      }
      #计算偏差
      CD_mat[i] = CD_function(U)
      WD_mat[i] = WD_function(U)
      MD_mat[i] = MD_function(U)
    }
    #结果1：输出结果矩阵，包含生成向量以及对应设计矩阵的CD,WD,MD
    result_mat = cbind(U_gene_mat,t(CD_mat),t(WD_mat),t(MD_mat))
    colnames(result_mat)[(s+1):(s+3)] = c('CD','WD','MD')
    #结果2：输出最小CD偏差对应的设计矩阵
    U_gene_CD = result_mat[which.min(result_mat[,s+1]),1:s]
    U_CD = matrix(0,nrow = n,ncol = s)
    U_CD[,1] = c(1:n)
    U_CD[1,]= U_gene_CD
    for(j in 2:s){
      U_CD[2:n,j] = mod(U_CD[2:n,1]*U_CD[1,j],n)
    }
    #结果3：输出最小WD偏差对应的设计矩阵
    U_gene_WD = result_mat[which.min(result_mat[,s+2]),1:s]
    U_WD = matrix(0,nrow = n,ncol = s)
    U_WD[,1] = c(1:n)
    U_WD[1,]= U_gene_WD
    for(j in 2:s){
      U_WD[2:n,j] = mod(U_WD[2:n,1]*U_WD[1,j],n)
    }
    #结果4：输出最小MD偏差对应的设计矩阵
    U_gene_MD = result_mat[which.min(result_mat[,s+3]),1:s]
    U_MD = matrix(0,nrow = n,ncol = s)
    U_MD[,1] = c(1:n)
    U_MD[1,]= U_gene_MD
    for(j in 2:s){
      U_MD[2:n,j] = mod(U_MD[2:n,1]*U_MD[1,j],n)
    }
    glpm_result = list(result_mat,U_CD,U_WD,U_MD)
    names(glpm_result) = c('result','U_CD','U_WD','U_MD')
    return(glpm_result)
  }else{
    print('m<s导致无法构造均匀设计！')
  }
}


###方幂好格子点法#######
##给定n和s####
n = 13;s = 3
fangmiglpm = function(n,s){
  library(schoolmath)
  #初始化均匀设计矩阵
  U = matrix(0,nrow = n,ncol = s)
  U[,1] = c(1:n)
  #构造备选正整数集A_ns
  #首先找到互质整数集##
  H_n = list()
  i = 1
  for(h in 1:n){
    if(gcd(n,h) == 1){
      H_n[i] = h
      i = i+1
    }
  }
  H_n = unlist(H_n)
  #从互质整数集中找到次数不小于s-1的正整数
  A_ns = list()
  i = 1
  for(h in H_n){
    U_gene_test = h^c(0:(s-1))
    if(length(unique(U_gene_test)) == length(U_gene_test)){
      A_ns[i] = h
      i = i+1
    }
  }
  A_ns = unlist(A_ns)
  #number表示一共由可以构造多少组均匀设计
  number = length(A_ns)
  #判断number是否大于0，否则无法构造设计矩阵
  if(number>=1){
    #初始化偏差矩阵,共可产生C^(m-1)_(s-1)个设计矩阵
    CD_mat = matrix(0,nrow = 1,ncol = number)
    WD_mat = matrix(0,nrow = 1,ncol = number)
    MD_mat = matrix(0,nrow = 1,ncol = number)
    #构造生成向量矩阵
    U_gene_mat = matrix(0,nrow = number,ncol = s)
    U_gene_mat[,1] = 1
    #使用事先定义的同余函数mod
    for(i in 1:number){
      U_gene_mat[i,2:s] =mod(A_ns[i]^c(1:(s-1)),n)
    }
    
    #构造设计矩阵
    for(i in 1:number){
      U[1,]= U_gene_mat[i,]
      #生成元a
      a = A_ns[i]
      for(j in 2:s){
        U[2:n,j] = mod(U[2:n,1]*a^(j-1),n)
      }
      #计算偏差
      CD_mat[i] = CD_function(U)
      WD_mat[i] = WD_function(U)
      MD_mat[i] = MD_function(U)
    }
    #结果1：输出结果矩阵，包含生成向量以及对应设计矩阵的CD,WD,MD
    result_mat = cbind(U_gene_mat,t(CD_mat),t(WD_mat),t(MD_mat))
    colnames(result_mat)[(s+1):(s+3)] = c('CD','WD','MD')
    #结果2：输出最小CD偏差对应的设计矩阵
    U_gene_CD = result_mat[which.min(result_mat[,s+1]),1:s]
    U_CD = matrix(0,nrow = n,ncol = s)
    U_CD[,1] = c(1:n)
    U_CD[1,]= U_gene_CD
    for(j in 2:s){
      U_CD[2:n,j] = mod(U_CD[2:n,1]*U_CD[1,j],n)
    }
    #结果3：输出最小WD偏差对应的设计矩阵
    U_gene_WD = result_mat[which.min(result_mat[,s+2]),1:s]
    U_WD = matrix(0,nrow = n,ncol = s)
    U_WD[,1] = c(1:n)
    U_WD[1,]= U_gene_WD
    for(j in 2:s){
      U_WD[2:n,j] = mod(U_WD[2:n,1]*U_WD[1,j],n)
    }
    #结果4：输出最小MD偏差对应的设计矩阵
    U_gene_MD = result_mat[which.min(result_mat[,s+3]),1:s]
    U_MD = matrix(0,nrow = n,ncol = s)
    U_MD[,1] = c(1:n)
    U_MD[1,]= U_gene_MD
    for(j in 2:s){
      U_MD[2:n,j] = mod(U_MD[2:n,1]*U_MD[1,j],n)
    }
    fangmiglpm_result = list(result_mat,U_CD,U_WD,U_MD)
    names(fangmiglpm_result) = c('result','U_CD','U_WD','U_MD')
    return(fangmiglpm_result)
  }else{
    print('number=0,无法构造设计矩阵！')
  }
}


#####切割法####
#需要通过好格子点法或者方幂好格子点法构造初始矩阵
#以MD为准则
#n为实际需要的设计矩阵行数
#n = 9
cutting = function(U_init,n){
  p = nrow(U_init);s = ncol(U_init)
  #构造不同准则下的偏差矩阵p*s规格
  CD_mat = matrix(0,nrow = p,ncol = s)
  WD_mat = matrix(0,nrow = p,ncol = s)
  MD_mat = matrix(0,nrow = p,ncol = s)
  #行排列
  for(l in 1:s){
    #记排列后的矩阵为U_init_l
    U_init_l = U_init[order(U_init[,l]),]
    #稍作变换，方便切割
    U_init_l_2 = rbind(U_init_l,U_init_l)
    #切割
    for(m in 1:p){
      #判断m与n的大小关系，作不同切割处理
      if(m<n){
        U_init_lm = U_init_l_2[c((p+m-n+1):(p+m)),]
      }else{
        U_init_lm = U_init_l[c((m-n+1):m),]
      }
      #切割后的矩阵根据排列顺序重记元素，记为U_lm
      U_lm = cbind(rank(U_init_lm[,1]),rank(U_init_lm[,2]),
                   rank(U_init_lm[,3]))
      CD_mat[m,l] = CD_function(U_lm)
      WD_mat[m,l] = WD_function(U_lm)
      MD_mat[m,l] = MD_function(U_lm)
    }
  }
  ##############################
  #找到最小CD对应设计
  index = which(CD_mat == min(CD_mat),arr.ind = T)
  if(dim(index)[1] == 1){
    l = index[2]
    m = index[1]
  }else{
    l = index[1,][2]
    m = index[1,][1]
  }
  #记排列后的矩阵为U_init_l
  U_init_l = U_init[order(U_init[,l]),]
  #稍作变换，方便切割
  U_init_l_2 = rbind(U_init_l,U_init_l)
  #切割
  #判断m与n的大小关系，作不同切割处理
  if(m<n){
    U_init_lm = U_init_l_2[c((p+m-n+1):(p+m)),]
  }else{
    U_init_lm = U_init_l[c((m-n+1):m),]
  }
  #切割后的矩阵根据排列顺序重记元素，记为U_lm
  U_CD = cbind(rank(U_init_lm[,1]),rank(U_init_lm[,2]),
               rank(U_init_lm[,3]))
  ###########################
  ##############################
  #找到最小WD对应设计
  index = which(WD_mat == min(WD_mat),arr.ind = T)
  if(dim(index)[1] == 1){
    l = index[2]
    m = index[1]
  }else{
    l = index[1,][2]
    m = index[1,][1]
  }
  #记排列后的矩阵为U_init_l
  U_init_l = U_init[order(U_init[,l]),]
  #稍作变换，方便切割
  U_init_l_2 = rbind(U_init_l,U_init_l)
  #切割
  #判断m与n的大小关系，作不同切割处理
  if(m<n){
    U_init_lm = U_init_l_2[c((p+m-n+1):(p+m)),]
  }else{
    U_init_lm = U_init_l[c((m-n+1):m),]
  }
  #切割后的矩阵根据排列顺序重记元素，记为U_lm
  U_WD = cbind(rank(U_init_lm[,1]),rank(U_init_lm[,2]),
               rank(U_init_lm[,3]))
  ###########################
  ##############################
  #找到最小MD对应设计
  index = which(MD_mat == min(MD_mat),arr.ind = T)
  if(dim(index)[1] == 1){
    l = index[2]
    m = index[1]
  }else{
    l = index[1,][2]
    m = index[1,][1]
  }
  #记排列后的矩阵为U_init_l
  U_init_l = U_init[order(U_init[,l]),]
  #稍作变换，方便切割
  U_init_l_2 = rbind(U_init_l,U_init_l)
  #切割
  #判断m与n的大小关系，作不同切割处理
  if(m<n){
    U_init_lm = U_init_l_2[c((p+m-n+1):(p+m)),]
  }else{
    U_init_lm = U_init_l[c((m-n+1):m),]
  }
  #切割后的矩阵根据排列顺序重记元素，记为U_lm
  U_MD = cbind(rank(U_init_lm[,1]),rank(U_init_lm[,2]),
               rank(U_init_lm[,3]))
  ###########################
  result = list(CD_mat,WD_mat,MD_mat,U_CD,U_WD,U_MD)
  names(result) = c('CD_mat','WD_mat','MD_mat',
                    'U_CD','U_WD','U_MD')
  return(result)
}

###调用函数
#例1
glpm_result = glpm(n=15,s=2)
#生成向量及对应偏差值
glpm_gene_devi = glpm_result$result
glpm_result[["U_CD"]]
glpm_result[["U_WD"]]
glpm_result[["U_MD"]]
#例2
#利用方幂glpm构造一个较大的均匀设计矩阵，选择依据为MD
U_init = fangmiglpm(17,3)
U_init$U_MD
#计算MD
MD_function(U_init$U_MD)
#切割法,n = 9
cutting_result = cutting(U_init$U_MD,9)
cutting_MD = cutting_result$MD_mat
cutting_result$U_MD
MD_function(cutting_result$U_MD)
