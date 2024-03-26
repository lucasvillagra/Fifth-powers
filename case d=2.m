function Hilbert_newforms(N)
        /*
        INPUT:
                - ``N`` : an integer

        OUTPUT:
                The space of classical newforms of level N as Hilbert newforms
        */

        QQ := RationalsAsNumberField();
        ZZ := Integers(QQ);
        M := HilbertCuspForms(QQ, N*ZZ);
        p := Min([p : p in PrimeFactors(N) | Valuation(N,p) eq 1]);
        A := QuaternionAlgebra(p*ZZ, InfinitePlaces(QQ) : Optimized);
        NewFs := NewSubspace(M : QuaternionOrder:=MaximalOrder(A));

        return  NewformDecomposition(NewFs);

end function;


//========================================================================================


function eliminate_newforms(k2 : Bound := 100, Newfs := [])
	/*
	INPUT:
		- ``kq`` : 1 or q
		- ``Bound`` : an upper bound of the primes we use in the elimination step
		- ``Newfs`` : a list of new forms

	OUTPUT:
		A list of primes for which the eliminations step fails.
	*/


	// We compute the newforms if they are not given


    if #Newfs eq 0 then
        if k2 eq 1 then
           	Newfs := Hilbert_newforms(2^8 * 3 * 5 * 11 * 17);
        elif k2 eq 2 then
            Newfs := Hilbert_newforms(3 * 5 * 11 * 17);
        end if;
	end if; 


	// We apply Proposition 2.1


	num := 0;	
	tHil := MakeType("ModFrmHil");
    Bfnews:=[];
	for fnew in Newfs do
   
			for q in PrimesInInterval(1,Bound) do
                if q notin [2,3,5,7,13] then
				
				    if Type(fnew) eq tHil then
                        ZZ := RingOfIntegers(BaseField(fnew));
					    aqfnew := HeckeEigenvalue(Eigenform(fnew),q*ZZ);
				    else
					    aqfnew := Coefficient(fnew[1],q);
				    end if;
				    prod := q * Norm((aqfnew^2 - (q+1)^2));

				    for r in [0..Floor(Sqrt(q))] do
					    prod *:= Norm((4*r^2 - aqfnew^2));
				    end for;

				    Append(~Bfnews,Integers()!prod);
                    g := GCD(Bfnews);
                    
                end if;		
			end for;


    g := GCD(Bfnews);
			if g ne 0 then
                print(PrimeFactors(g));
			else
				print "We have a zero";
			end if;
end for;
return 0;

end function;
