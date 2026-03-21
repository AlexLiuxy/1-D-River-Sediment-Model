function res = Fe_bc(Fea,Feb)
global Feinit
  res = [ Fea(1)-Feinit 
          Feb(2) ];
 
end