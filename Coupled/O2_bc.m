
function res = O2_bc(O2a,O2b)
global O2init

  res = [ O2a(1)-O2init 
          O2b(2) ];
 
end