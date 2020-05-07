%%
% mdfData.Time
figure()
title('experiment v1.1')
subplot(2,2,1)
plot(mdfData.Time,mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightL_mm__In1)
title('T1')

subplot(2,2,2)
plot(mdfData.Time,mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightM_mm__In1)
title('T3')

subplot(2,2,3)
plot(mdfData.Time,mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightR_mm__In1)
title('T2')