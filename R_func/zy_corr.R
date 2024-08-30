# 代码虽然罗嗦，但是不想改 
zy_corr = function(x, y=NA, method="spearman", type="matrix", step=100, na_rm=T){
  # 计算行与行之间的相关系数
  # type -> 返回matrix或长列表 c("matrix", "long")
  # 长列表只返回a,b, 不返回b,a的相关性
  # step 表示每隔多少个回显一条信息
  num_x =  nrow(x)
  name_x =  rownames(x)
  # y 不为空
  if(!is.na(y)){
    name_y = rownames(y)
    num_y = nrow(y)
    # 如果是matrix
    if(type=="matrix"){
      p_result = cor_result = matrix(NA,nrow=num_x, ncol=num_y, dimnames=list(name_x, name_y))
      for(i in 1:num_x){
        if( i %% step == 0){message("num: ",i)}
        name_i = name_x[i]
        message(name_i)
        for(j in 1:num_y){
          name_j = name_y[j]
          ct = cor.test(as.numeric(x[name_i,]), as.numeric(y[name_j,]), emthod=method, na.rm=na_rm)
          cor_result[name_i,name_j] = ct$estimate
          p_result[name_i,name_j] = ct$p.value
        }
      }
    }else if(type=="long"){
      c = 1
      result = matrix(NA,nrow=num_x*num_y, ncol=4, 
                      dimnames=list(NULL,c("name_a","name_b","corr","pval")))
      for(i in 1:num_x){
        if( i %% step == 0){message("num: ",i)}
        name_i = name_x[i]
        for(j in 1:num_y){
          name_j = name_y[j]
          ct = cor.test(as.numeric(x[name_i,]), as.numeric(y[name_j,]), method=method, na.rm=na_rm)
          result[c,]  = c(name_i, name_j, ct$estimate, ct$p.value)
          c = c+1
        }
      }
    }
  }else{# 只有x
    if(type=="matrix"){
      p_result = cor_result = matrix(NA,nrow=num_x, ncol=num_x,
                      dimnames=list(name_x, name_x))
      for(i in 1:num_x){
        if( i %% step == 0){message("num: ",i)}
        name_i = name_x[i]
        for(j in i:num_x){
          name_j = name_x[j]
          ct = cor.test(as.numeric(x[name_i,]), as.numeric(x[name_j,]), method=method, na.rm=na_rm)
          cor_result[name_i,name_j] = cor_result[name_j,name_i] = ct$estimate
          p_result[name_i,name_j] = p_result[name_j, name_i]= ct$p.value
        }
      }
    }else if(type=="long"){
      result = matrix(NA,nrow=num_x*(num_x+1)/2, ncol=4,
                      dimnames=list(NULL,c("name_a","name_b","corr","pval")))
      c = 1
      for(i in 1:num_x){
        if( i %% step == 0){message("num: ",i)}
        name_i = name_x[i]
        message(name_i)
        for(j in i:num_x){
          name_j = name_x[j]
          ct = cor.test(as.numeric(x[name_i,]), as.numeric(x[name_j,]), method=method, na.rm=na_rm)
          result[c,]  = c(name_i, name_j, ct$estimate, ct$p.value)
          c = c+1
        }
      }
    }
  }
  if(type=="long"){
    return(result)
  }else{return(list(p=p_result, c=cor_result))}
}

zy_corr_row = zy_corr
