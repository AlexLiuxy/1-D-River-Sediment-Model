
function res = CH4_bc(CH4a,CH4b)
global CH4init
  res = [ CH4a(1)-CH4init 
          CH4b(2) ];
 
end