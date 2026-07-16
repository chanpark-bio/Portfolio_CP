#!/bin/bash

echo "================================================="
echo "  CP In Silico Hub - New Project Wizard"
echo "================================================="

if [ "$#" -eq 3 ]; then
    PROJECT_NAME=$1
    TARGET_NAME=$2
    PDB_ID=$3
else
    echo "새로운 프로젝트 생성을 위한 정보를 입력해 주십시오."
    echo ""
    read -p "1. 프로젝트명을 입력하세요 (예: Cancer_DeNovo) : " PROJECT_NAME
    read -p "2. 타겟 단백질명을 입력하세요 (예: PDL1)         : " TARGET_NAME
    read -p "3. 타겟 PDB ID를 입력하세요 (예: Q9NZQ7)         : " PDB_ID
    echo ""
fi

echo "⚙️  [$PROJECT_NAME\_$TARGET_NAME] 프로젝트 인프라를 구축합니다..."
python3 01_Library/Scripts/DeNovo_Pipeline/00_init_workspace.py \
    --project "$PROJECT_NAME" \
    --target "$TARGET_NAME" \
    --pdb_id "$PDB_ID"
