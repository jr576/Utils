package io.jr576;

import java.util.ArrayList;
import java.util.List;

public class BinaryIndexedTree {
    final long[] binaryIndexedTree;
    final int upperbound;

    public BinaryIndexedTree(int upperbound) {
        this.upperbound = upperbound;
        binaryIndexedTree = new long[upperbound + 1];
    }

    public long getSum(int index) {
        long sum = 0;
        while (index > 0) {
            sum += binaryIndexedTree[index];
            index -= Integer.lowestOneBit(index);
        }
        return sum;
    }

    public long getSumMod(int index, long mod){
        long sum = 0;
        while (index > 0) {
            sum = (sum + binaryIndexedTree[index]) % mod;
            index -= Integer.lowestOneBit(index);
        }
        return sum;
    }

    public void add(int index, long val) {
        while (index <= upperbound) {
            binaryIndexedTree[index] += val;
            index += Integer.lowestOneBit(index);
        }
    }

    public void addMod(int index, long val, long mod) {
        while (index <= upperbound) {
            binaryIndexedTree[index] = (binaryIndexedTree[index]+val) % mod;
            index += Integer.lowestOneBit(index);
        }
    }

    @Override
    public String toString() {
        List<Long> tempList = new ArrayList<>();
        for (int i = 0; i <= upperbound; i++)
            tempList.add(getSum(i));
        return tempList.toString();
    }
}

