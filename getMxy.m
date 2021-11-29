function Mxy = getMxy(M)
    
    Mx = squeeze(M(:,1,:));
    My = squeeze(M(:,2,:));
    Mxy = Mx + 1i*My;

end

