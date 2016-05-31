ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end 