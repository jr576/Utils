package io.jr576.utils;


import java.util.Arrays;
import java.util.Objects;
import java.util.function.LongBinaryOperator;
import java.util.function.LongUnaryOperator;
import java.util.stream.IntStream;

public class Dirichlet {
    public long[] values;
    public final long[] quotients;
    public final long upperbound;
    public final int upperboundSqrt;
    final boolean removeOne;
    final long breakPoint; //upperbound^(2/3)

    public Dirichlet(long upperbound) {
        this.upperbound = upperbound;
        this.upperboundSqrt = (int) Math.sqrt(upperbound);
        this.removeOne = upperbound / upperboundSqrt == upperboundSqrt;
        values = new long[2 * upperboundSqrt + (this.removeOne ? 0 : 1)];
        quotients = new long[2 * upperboundSqrt + (this.removeOne ? 0 : 1)];
        for (int d = 1; d <= upperboundSqrt; d++)
            quotients[d] = d;
        for (int d = upperboundSqrt; d >= 1; d--)
            quotients[2 * upperboundSqrt + (this.removeOne ? 0 : 1) - d] = upperbound / d;
        this.breakPoint = largestQuotientAtMost((long) Math.pow(upperbound, 2.0/3));
    }

    public Dirichlet(long upperbound, LongUnaryOperator function) {
        this(upperbound);
        for (int index = 1; index < quotients.length; index++) {
            long quotient = quotients[index];
            values[index] = function.applyAsLong(quotient);
        }
    }

    public Dirichlet clone() {
        return new Dirichlet(this);
    }

    private Dirichlet(Dirichlet copyFieldsOf) {
        this.upperbound = copyFieldsOf.upperbound;
        this.upperboundSqrt = copyFieldsOf.upperboundSqrt;
        this.removeOne = copyFieldsOf.removeOne;
        this.values = new long[copyFieldsOf.values.length];
        this.quotients = copyFieldsOf.quotients.clone();
        this.breakPoint = copyFieldsOf.breakPoint;
    }

    public int quotientToIndex(long quotient) {
        if (quotient <= upperboundSqrt)
            return (int) quotient;
        else
            return (int) (2 * upperboundSqrt - upperbound / quotient) + (removeOne ? 0 : 1);
    }

    public long largestQuotientAtMost(long integer) {
        return quotients[quotientToIndex(integer)];
    }

    public long get(long quotient) {
        return values[quotientToIndex(quotient)];
    }

    public void put(long quotient, long value) {
        values[quotientToIndex(quotient)] = value;
    }

    public void put(long first, long last, long value) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] = value;
            index--;
            quotient = quotients[index];
        }
    }

    public void put(long first, long last, LongUnaryOperator function) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] = function.applyAsLong(quotient);
            index--;
            quotient = quotients[index];
        }
    }

    public void add(long quotient, long valueToAdd) {
        values[quotientToIndex(quotient)] += valueToAdd;
    }

    public void add(long first, long last, long valueToAdd) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] += valueToAdd;
            index--;
            quotient = quotients[index];
        }
    }

    public void add(long first, long last, LongUnaryOperator toAddFunction) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] += toAddFunction.applyAsLong(quotient);
            index--;
            quotient = quotients[index];
        }
    }

    public void addMod(long quotient, long valueToAdd, long mod) {
        int quotientIndex = quotientToIndex(quotient);
        values[quotientIndex] = (values[quotientIndex] + valueToAdd) % mod;
    }

    public void addMod(long first, long last, long valueToAdd, long mod) {
        long quotient = last;
        int quotientIndex = quotientToIndex(last);
        while (quotient >= first) {
            values[quotientIndex] = (valueToAdd + values[quotientIndex]) % mod;
            quotientIndex--;
            quotient = quotients[quotientIndex];
        }
    }

    public void addMod(long first, long last, LongUnaryOperator toAddFunction, long mod) {
        long quotient = last;
        int quotientIndex = quotientToIndex(last);
        while (quotient >= first) {
            values[quotientIndex] = (values[quotientIndex] + toAddFunction.applyAsLong(quotient)) % mod;
            quotientIndex--;
            quotient = quotients[quotientIndex];
        }
    }

    public void summatory() {
        for (int quotientIndex = 1; quotientIndex < quotients.length; quotientIndex++)
            values[quotientIndex] += values[quotientIndex - 1];
    }

    public Dirichlet square() {
        Dirichlet square = new Dirichlet(this);
        IntStream.iterate(0, i -> i < quotients.length, i -> i + 1)
                 .parallel()
                 .forEach(i -> squareHelper(quotients[i], square));
        return square;
    }

    private Void squareHelper(long quotient, Dirichlet square) {
        long quotientSqrt = (long) Math.sqrt(quotient);
        long value = -this.get(quotientSqrt) * this.get(quotientSqrt);
        for (int d = 1; d <= quotientSqrt; d++)
            value += 2 * (this.get(d) - this.get(d - 1)) * this.get(quotient / d);
        square.put(quotient, value);
        return null;
    }

    public Dirichlet squareMod(long mod) {
        Dirichlet square = new Dirichlet(this);
        IntStream.iterate(0, i -> i < quotients.length, i -> i + 1)
                 .parallel()
                 .forEach(i -> squareModHelper(quotients[i], square, mod));
        return square;
    }

    private void squareModHelper(long quotient, Dirichlet square, long mod) {
        long quotientSqrt = (long) Math.sqrt(quotient);
        long value = -this.get(quotientSqrt) * this.get(quotientSqrt) % mod;
        for (int d = 1; d <= quotientSqrt; d++)
            value = (value + 2 * (this.get(d) - this.get(d - 1)) * this.get(quotient / d)) % mod;
        square.put(quotient, Math.floorMod(value, mod));
    }

    public Dirichlet multiply(Dirichlet multiplyWith) {
        Dirichlet product = new Dirichlet(this);
        IntStream.iterate(0, i -> i < quotients.length, i -> i + 1)
                 .parallel()
                 .forEach(i -> multiplyHelper(quotients[i], product, multiplyWith));
        return product;
    }

    public Dirichlet multiplyTwo(Dirichlet multiplyWith) {
        Dirichlet product = new Dirichlet(this);
        long breakPoint = largestQuotientAtMost((long) Math.pow(upperbound, 2.0/3));
        for (int i = 1; quotients[i] <= breakPoint; i++) {
            for (int j = 1; quotients[j] <= breakPoint / quotients[i]; j++)
                product.add(quotients[i] * quotients[j], (this.values[i] - this.values[i - 1]) * (multiplyWith.values[j] - multiplyWith.values[j - 1]));
            product.values[i] = product.values[i] + product.values[i - 1];
        }
        for (int i = quotients.length - 1; quotients[i] > breakPoint; i--)
            multiplyHelper(quotients[i], product, multiplyWith);
        return product;
    }

    public Dirichlet multiplyModTwo(Dirichlet multiplyWith, long mod) {
        Dirichlet product = new Dirichlet(this);
        long breakPoint = largestQuotientAtMost((long) Math.pow(upperbound, 2.0/3));
        for (int i = 1; quotients[i] <= breakPoint; i++) {
            for (int j = 1; quotients[j] <= breakPoint / quotients[i]; j++)
                product.addMod(quotients[i] * quotients[j], (this.values[i] - this.values[i - 1]) * (multiplyWith.values[j] - multiplyWith.values[j - 1]), mod);
            product.values[i] = product.values[i] + product.values[i - 1];
        }
        for (int i = quotients.length - 1; quotients[i] > breakPoint; i--)
            multiplyModHelper(quotients[i], product, multiplyWith, mod);
        return product;
    }

    private void multiplyHelper(long quotient, Dirichlet product, Dirichlet multiplyWith) {
        long quotientSqrt = (long) Math.sqrt(quotient);
        long value = -this.get(quotientSqrt) * multiplyWith.get(quotientSqrt);
        for (int i = 1; i <= quotientSqrt; i++)
            value += (this.get(i) - this.get(i - 1)) * multiplyWith.get(quotient / i)
                     + this.get(quotient / i) * (multiplyWith.get(i) - multiplyWith.get(i - 1));
        product.put(quotient, value);
    }

    public Dirichlet multiplyMod(Dirichlet multiplyWith, long mod) {
        Dirichlet product = new Dirichlet(this);
        IntStream.iterate(0, i -> i < quotients.length, i -> i + 1)
                 .parallel()
                 .forEach(i -> multiplyModHelper(quotients[i], product, multiplyWith, mod));
        return product;
    }

    private void multiplyModHelper(long quotient, Dirichlet multiplyWith, Dirichlet product, long mod) {
        long quotientSqrt = (long) Math.sqrt(quotient);
        long value = -this.get(quotientSqrt) * multiplyWith.get(quotientSqrt) % mod;
        for (int i = 1; i <= quotientSqrt; i++)
            value = (value + (this.get(i) - this.get(i - 1)) * multiplyWith.get(quotient / i)
                     + this.get(quotient / i) * (multiplyWith.get(i) - multiplyWith.get(i - 1))) % mod;
        product.put(quotient, Math.floorMod(value, mod));
    }

    public static Dirichlet one(long upperbound) {
        return new Dirichlet(upperbound, i -> 1);
    }

    public Dirichlet pow(long exponent) {
        if (exponent < 1) {
            if (exponent < 0)
                return null;
            return Dirichlet.one(upperbound);
        }
        Dirichlet base = this;
        while (exponent % 2 == 0) {
            exponent /= 2;
            base = base.square();
        }
        Dirichlet result = base;
        exponent /= 2;
        while (exponent > 0) {
            base = base.square();
            if (exponent % 2 == 1)
                result = result.multiply(base);
            exponent /= 2;
        }
        return result;
    }

    public Dirichlet powMod(long exponent, long mod) {
        if (exponent < 1) {
            if (exponent < 0)
                return null;
            return Dirichlet.one(upperbound);
        }
        Dirichlet base = this;
        while (exponent % 2 == 0) {
            exponent /= 2;
            base = base.squareMod(mod);
        }
        Dirichlet result = base;
        exponent /= 2;
        while (exponent > 0) {
            base = base.squareMod(mod);
            if (exponent % 2 == 1)
                result = result.multiplyMod(base, mod);
            exponent /= 2;
        }
        return result;
    }

    public void map(LongUnaryOperator function) {
        for (int index = values.length - 1; index >= 1; index--)
            values[index] = function.applyAsLong(values[index]);
    }

    public void map(LongBinaryOperator function) {
        for (int index = values.length - 1; index >= 1; index--)
            values[index] = function.applyAsLong(quotients[index], values[index]);
    }

    public void map(long first, long last, LongUnaryOperator function) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] = function.applyAsLong(values[index]);
            index--;
            quotient = quotients[index];
        }
    }

    public void map(long first, long last, LongBinaryOperator function) {
        long quotient = last;
        int index = quotientToIndex(last);
        while (quotient >= first) {
            values[index] = function.applyAsLong(quotients[index], values[index]);
            index--;
            quotient = quotients[index];
        }
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (int quotientIndex = 0; quotientIndex <= this.quotients.length - 1; quotientIndex++)
            stringBuilder.append(this.quotients[quotientIndex]).append("=").append(this.values[quotientIndex]).append(", ");
        stringBuilder.delete(stringBuilder.length() - 2, stringBuilder.length());
        return "{" + stringBuilder + "}";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Dirichlet dirichlet = (Dirichlet) o;
        return upperbound == dirichlet.upperbound && upperboundSqrt == dirichlet.upperboundSqrt && removeOne == dirichlet.removeOne && Arrays.equals(values, dirichlet.values) && Arrays.equals(quotients, dirichlet.quotients);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(upperbound, upperboundSqrt, removeOne);
        result = 31 * result + Arrays.hashCode(values);
        result = 31 * result + Arrays.hashCode(quotients);
        return result;
    }


}

