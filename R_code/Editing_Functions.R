#Filtering criteria
A/T reference, in REDIportal, > 20 reads, <.99 editing, not .5 editing
#Select ADAR editing sites
filtered <- filter(reads, REF == "T" | REF == "A")
#attach editing frequency column
filtered <- mutate(filtered, Editing.rate = 0)
#Calculate editing frequency for each site: distinguishes between A->G and T->c
for (n in 1:nrow(filtered)){
  if (filtered$REF[n] == "A"){
    filtered$Editing.rate[n] = filtered$G.count[n]/filtered$TOTAL[n]
  }
  else{
    filtered$Editing.rate[n] = filtered$C.count[n]/filtered$TOTAL[n]
  }
}


return(filtered)
}


}

#Function for attaching mean editing rate assuming that the data frame contains position (POS) column

left_join(Sampledataframe, 
          Sampledataframe %>% group_by(POS) %>% 
            summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
            ungroup()



















