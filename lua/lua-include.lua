util.biot = util.biot or {}


function util.biot.PrintDebug(text)
  print("LUA/Poroelasticity:\t"..text)
end


function util.biot.CheckAssertions()
  if not util.json then util.biot.PrintDebug("FATAL: util.json not found"); quit() end
end

-- TODO: Go C++
function util.biot.CreateApproxSpace(dom, dim, uorder, porder)
print("Create ApproximationSpace... ")
local approxSpace = ApproximationSpace(dom) 
approxSpace:add_fct("p", "Lagrange", porder) 

if false then 
  -- Does not work due to registration issues in SmallStrain mechanics.
  uorder=1
  approxSpace:add_fct("ux", "mini", 1)          
  approxSpace:add_fct("uy", "mini", 1)  
else
  local utype = "Lagrange" --"Lagrange"  --"mini" -- "Lagrange"
  approxSpace:add_fct("ux", utype, uorder)          
  approxSpace:add_fct("uy", utype, uorder)   
  if (dim==3) then approxSpace:add_fct("uz", utype, uorder) end
end

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_layout_statistic()
approxSpace:print_local_dof_statistic(2)              
print("... done!")
return approxSpace
end 


function util.biot.CreateDefaultErrorEst(dim)
  local biotErrorEst 

  biotErrorEst = ScaledGridFunctionEstimator()
  biotErrorEst:add(L2ComponentSpace("p", 2))        -- L2 norm for p, 2nd order quadrature

  biotErrorEst:add(H1SemiComponentSpace("ux", 4))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4)) 
  if (dim==3) then biotErrorEst:add(H1SemiComponentSpace("uz", 4)) end
  return biotErrorEst
end

