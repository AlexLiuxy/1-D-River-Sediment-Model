
function res = SO4_bc(SO4a,SO4b)
global SO4init
  res = [ SO4a(1)-SO4init 
          SO4b(2) ];
 
end