function y = gen_no(b,B) % given branch b from B, output generation no.
    if b(1,1)==B(1,1) % root, define as generation 1
        y=1;
    else % recursion using parent branch
        b=B(B(:,2)==b(1,1),:);
        y=1+gen_no(b,B);
    end
end