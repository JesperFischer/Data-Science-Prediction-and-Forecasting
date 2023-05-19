#utility functions



get_3darray = function(data,variable){
  if(variable == "percept"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    
    return(qq2)
  }else if(variable == "percept_bin"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    return(qq2)
  }else if(variable == "pred"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    return(qq2)
  }else if(variable == "stim"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    return(qq2)
    
  }else if(variable == "cue"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    return(qq2)
    
  }else if(variable == "u"){
    qq  = as.matrix(data %>% filter(session == 1) %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
    qq1  = as.matrix(data %>% filter(session == 2) %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
    qq2 = simplify2array(list(qq,qq1))
    return(qq2)
  }else{
    print("error")
  }
}
