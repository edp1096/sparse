## Memo
* Singular : 역행렬이 없는, 행렬식이 0인 행렬
* Singleton : 행 또는 열에 단 하나의 비영(non-zero)원소만 존재
* Fill-in : LU 또는 Cholesky 등의 분해 중, 본래 0인 위치에 새로운 비영원소가 생기는 것
* Pivot : 분해시 현재 선택/처리되는 대각원소
* Partition : 큰 행렬을 더 작은 부분 행렬들로 나누는 것. 블록 행렬 연산이나 병렬 처리 목적
* TrashCan : 희소 행렬 연산 중에 필요 없는 원소들의 임시로 저장하는 공간. 메모리 관리와 가비지 컬렉션 목적

## TODO
* app
    * mat0, matrix4000 create with size 0 -> Factor -> panic
