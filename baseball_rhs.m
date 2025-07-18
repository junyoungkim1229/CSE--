% 야구공의 운동 방정식 정의
function dsdt = baseball_rhs(s, F, B, g, phi, omega)
    % 상태 벡터: s = [x; y; z; vx; vy; vz]
    x = s(1); y = s(2); z = s(3);
    vx = s(4); vy = s(5); vz = s(6);
    
    V = sqrt(vx^2 + vy^2 + vz^2);
    f = F(V);

    dx = vx;
    dy = vy;
    dz = vz;

    dvx = -f*V*vx + B*omega*(vz*sin(phi) - vy*cos(phi));
    dvy = -f*V*vy + B*omega*vx*cos(phi);
    dvz = -g - f*V*vz - B*omega*vx*sin(phi);

    dsdt = [dx; dy; dz; dvx; dvy; dvz];
end
