void merge
(
    DynamicList<label>& list,
    const label low,
    const label mid,
    const label high,
    const volScalarField& solT
)
{
    label size1 = mid - low + 1;
    label size2 = high - mid;

    label sub1[size1], sub2[size2];
    scalar subSol1[size1], subSol2[size2];

    for(label i=0; i<size1; i++)
    {
        sub1[i]=list[low+i];
        label cellID = list[low+i];

        subSol1[i]=solT[cellID];
    }

    for(label j=0; j<size2; j++)
    {
        sub2[j]=list[mid+1+j];
        label cellID = list[mid+1+j];
      
        subSol2[j]=solT[cellID];
    }

    label i = 0;
    label j = 0;
    label k = low;

    while(i<size1 && j<size2)
    {
        if(subSol1[i] <= subSol2[j])
        {
            list[k] = sub1[i];
            i++;
        }
        else
        {
            list[k] = sub2[j];
            j++;
        }
        k++;
    }

    while (i < size1)
    {
        list[k] = sub1[i];
        i++;
        k++;
    }

    while (j < size2)
    {
        list[k] = sub2[j];
        j++;
        k++;
    }
}

void mergeSort
(
    DynamicList<label>& list,
    const label low,
    const label high,
    const volScalarField& solT
)
{ 
    if (low<high)
    {
        label mid = (low + high - 1)/2;
        mergeSort(list, low, mid, solT);
        mergeSort(list, mid+1, high, solT);

        merge(list, low, mid, high, solT);
    }
}

int sortDynamicList(DynamicList<label>& list, const volScalarField& solT)
{
    label size = list.size();
    mergeSort(list, 0, size - 1, solT);

    return 0;
} 
