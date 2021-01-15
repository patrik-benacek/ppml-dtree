get_crch_model <- function(train, target){
    if (target == 't2m'){ 
        model <- crch(obs~fc_mean | I(fc_std^2), data = train, dist="gaussian"  , link.scale = "quad", type='crps')
    }
    if (target == 'prec24'){ 
        model <- crch(obs~fc_mean | log(fc_std), data = train, dist = 'logistic', link.scale = "log" , left = 0, type='crps') # crch has better scores than trch
        #model <- crch(obs~fc_mean | log(fc_std), data = train, dist = 'logistic', link.scale = "log" , left = 0, truncated = TRUE, subset = obs>0) 
    }
    return(model)
}

get_crps_scores <- function(eval, target){
    if (target == 't2m'){ 
        score_raw = crps(eval$obs, location = eval$fc_mean, scale = eval$fc_std, family='normal')
        score_ngr = crps(eval$obs, location = eval$pp_mean, scale = eval$pp_std, family='normal')
    }
    if (target == 'prec24'){ 
        score_raw = crps(eval$obs, location = eval$fc_mean, scale = eval$fc_std, lower=0, upper=Inf, family='clogis')
        score_ngr = crps(eval$obs, location = eval$pp_mean, scale = eval$pp_std, lower=0, upper=Inf, family='clogis')
    }
    return(list(score_raw, score_ngr))
}
