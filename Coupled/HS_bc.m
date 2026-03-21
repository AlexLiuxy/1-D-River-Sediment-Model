
function res = HS_bc(HSa,HSb)
global HSinit
  res = [ HSa(1)-HSinit 
          HSb(2) ];
 
end