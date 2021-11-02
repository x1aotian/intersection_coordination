function plot_1(T,V,Nv)
    dset = ['--or';'--og';'--ob';'--oc';'--om';'--oy';'--ok';];
    for i = 1:Nv
        plot(T,V(i,:),dset(i));
        hold on;
    end
end