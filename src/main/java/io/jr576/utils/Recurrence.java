package io.jr576.utils;

//initial should consist of the initial values [a_0, a_1,...]
//recurrence array should be of form [c_{n-1}, c_{n-2},...] for a recurrence a_n = c_{n-1}a_{n-1}+c_{n-2}a_{n-2}+...
public class Recurrence {
    private final int degree;
    private final long[] initialValues;
    private final long[] recurrence;

    public Recurrence(long[] initial, long[] recurrence) {
        if(initial.length < recurrence.length)
            throw new IllegalArgumentException("InitialValues need to be at least as long as recurrence array");
        if(initial.length > recurrence.length) {
            this.initialValues = new long[recurrence.length];
            for(int i = 0; i < recurrence.length; i++)
                this.initialValues[i] = initial[i];
        } else this.initialValues = initial;
        this.recurrence = recurrence;
        this.degree = recurrence.length;
    }

    public long getTerm(long n) {
        int baseStart = 1;
        int polStart = 0;
        while(2*baseStart < degree && n > 0) {
            if(n%2 == 1)
                polStart += baseStart;
            baseStart *= 2;
            n /= 2;
        }
        long[] pol = new long[2*degree - 1];
        long[] base = new long[2*degree - 1];
        pol[polStart] = 1;
        base[baseStart] = 1;
        if(baseStart >= degree) reduce(base);
        if(polStart >= degree) reduce(pol);
        while(n > 0) {
            if(n%2 == 1)
                pol = multiply(pol, base);
            base = multiply(base, base);
            n /= 2;
        }
        long term = 0;
        for(int i = 0; i < degree; i++)
            term = (term + initialValues[i]*pol[i]);
        return term;
    }

    public long getTerm(long n, long mod) {
        if(n < initialValues.length)
            return initialValues[(int)n];
        int baseStart = 1;
        int polStart = 0;
        while(2*baseStart < degree && n > 0) {
            if(n%2 == 1)
                polStart += baseStart;
            baseStart *= 2;
            n /= 2;
        }
        long[] pol = new long[2*degree - 1];
        long[] base = new long[2*degree - 1];
        pol[polStart] = 1;
        base[baseStart] = 1;
        if(baseStart >= degree) reduce(base);
        if(polStart >= degree) reduce(pol);
        while(n > 0) {
            if(n%2 == 1)
                pol = multiply(pol, base, mod);
            base = multiply(base, base, mod);
            n /= 2;
        }
        long term = 0;
        for(int i = 0; i < degree; i++)
            term = (term + initialValues[i]*pol[i])%mod;
        return term;
    }

    private long[] multiply(long[] polA, long[] polB) {
        long[] product = new long[2*degree - 1];
        for(int x = 0; x < degree; x++)
            if(polA[x] != 0)
                for(int y = 0; y < degree; y++)
                    product[x + y] = (product[x + y] + polA[x]*polB[y]);
        reduce(product);
        return product;
    }

    private long[] multiply(long[] polA, long[] polB, long mod) {
        long[] product = new long[2*degree - 1];
        for(int x = 0; x < degree; x++)
            if(polA[x] != 0)
                for(int y = 0; y < degree; y++)
                    product[x + y] = (product[x + y] + polA[x]*polB[y])%mod;
        reduce(product, mod);
        return product;
    }
    
    private void reduce(long[] pol) { //reduction modulo characteristic polynomial
        for(int i = 2*degree - 2; i >= degree; i--) {
            if(pol[i] != 0)
                for(int j = 1; j <= degree; j++)
                    pol[i - j] = (pol[i - j] + recurrence[j - 1]*pol[i]);
            pol[i] = 0;
        }
    }

    private void reduce(long[] pol, long mod) { //reduction modulo characteristic polynomial with mod
        for(int i = 2*degree - 2; i >= degree; i--) {
            if(pol[i] != 0)
                for(int j = 1; j <= degree; j++)
                    pol[i - j] = (pol[i - j] + recurrence[j - 1]*pol[i])%mod;
            pol[i] = 0;
        }
    }
}

