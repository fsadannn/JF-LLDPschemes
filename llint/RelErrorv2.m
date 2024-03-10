function E = RelErrorv2(Y,Z,n_parts)

% Y : vector con valores exactos
% Z : vector con valores approximados

dim=size(Y,1).*size(Y,2);
y=reshape(Y,dim,1);
z=reshape(Z,dim,1);
index=find(abs(y)>eps);
% E = max(max(abs((y(index)-z(index))./y(index))));
y = y(index);
z= z(index);
E = 0;
idd = linspace(1,length(y),n_parts);
for i=2:length(idd)
    idx = floor(idd(i));
    idx_old = floor(idd(i-1));
   
    E = max( max(max(abs((y(idx_old:idx)-z(idx_old:idx))./y(idx_old:idx)))),E);
end

end