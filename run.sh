#!/bin/bash

source 00_System/activate_engines.sh

echo "================================================="
echo "  CP In Silico Hub - Interactive Auto Launcher"
echo "================================================="

read -p "▶ 프로젝트명을 입력하세요 (예: Obesity_Denovo): " PROJECT_NAME
read -p "▶ 타겟명을 입력하세요 (예: PDL1): " TARGET_NAME

WORKSPACE_DIR="03_Workspace/${PROJECT_NAME}_${TARGET_NAME}"
CONFIG_PATH="${WORKSPACE_DIR}/00_Master_Control/config.yaml"

if [ ! -f "$CONFIG_PATH" ]; then
    echo " [FATAL ERROR] 설정 파일이 없습니다. 오타를 확인하거나 init.sh를 실행하세요: $CONFIG_PATH"
    exit 1
fi

python3 00_System/interactive_wizard.py "$PROJECT_NAME" "$TARGET_NAME"
