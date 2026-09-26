// Trajectory Comparison Web Interface

class ComparisonManager {
    constructor() {
        this.trajectories = [];
        this.initializeElements();
        this.attachEventListeners();
    }

    initializeElements() {
        this.uploadZone = document.getElementById('uploadZone');
        this.fileInput = document.getElementById('fileInput');
        this.trajectoryList = document.getElementById('trajectoryList');
        this.trajCount = document.getElementById('trajCount');
        this.compareBtn = document.getElementById('compareBtn');
        this.resultsArea = document.getElementById('resultsArea');
        this.analysisType = document.getElementById('analysisType');
        this.selection = document.getElementById('selection');
        this.showStats = document.getElementById('showStats');
        this.alignTrajectories = document.getElementById('alignTrajectories');
    }

    attachEventListeners() {
        // Upload zone events
        this.uploadZone.addEventListener('click', () => this.fileInput.click());
        this.fileInput.addEventListener('change', (e) => this.handleFileSelect(e));

        // Drag and drop
        this.uploadZone.addEventListener('dragover', (e) => {
            e.preventDefault();
            this.uploadZone.classList.add('dragover');
        });

        this.uploadZone.addEventListener('dragleave', () => {
            this.uploadZone.classList.remove('dragover');
        });

        this.uploadZone.addEventListener('drop', (e) => {
            e.preventDefault();
            this.uploadZone.classList.remove('dragover');
            this.handleFileSelect({ target: { files: e.dataTransfer.files } });
        });

        // Compare button
        this.compareBtn.addEventListener('click', () => this.runComparison());
    }

    async handleFileSelect(event) {
        const files = Array.from(event.target.files);

        for (const file of files) {
            await this.uploadTrajectory(file);
        }

        this.fileInput.value = '';
    }

    async uploadTrajectory(file) {
        const formData = new FormData();
        formData.append('file', file);

        try {
            const response = await fetch('/api/comparison/upload', {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                throw new Error(`Upload failed: ${response.statusText}`);
            }

            const data = await response.json();
            this.addTrajectory({
                id: data.file_id,
                name: file.name,
                path: data.path,
                size: this.formatFileSize(file.size)
            });

            this.showMessage('success', `Uploaded: ${file.name}`);
        } catch (error) {
            this.showMessage('error', `Failed to upload ${file.name}: ${error.message}`);
        }
    }

    addTrajectory(trajectory) {
        this.trajectories.push(trajectory);
        this.renderTrajectoryList();
        this.updateCompareButton();
    }

    removeTrajectory(id) {
        this.trajectories = this.trajectories.filter(t => t.id !== id);
        this.renderTrajectoryList();
        this.updateCompareButton();

        // Delete from server
        fetch(`/api/comparison/upload/${id}`, { method: 'DELETE' })
            .catch(err => console.error('Failed to delete file:', err));
    }

    renderTrajectoryList() {
        this.trajCount.textContent = this.trajectories.length;

        if (this.trajectories.length === 0) {
            this.trajectoryList.innerHTML = '<p style="color: #999; text-align: center; padding: 20px;">No trajectories uploaded</p>';
            return;
        }

        this.trajectoryList.innerHTML = this.trajectories.map(traj => `
            <div class="trajectory-item" data-id="${traj.id}">
                <button class="remove-btn" onclick="comparisonManager.removeTrajectory('${traj.id}')">×</button>
                <div class="trajectory-name">${traj.name}</div>
                <div class="trajectory-info">Size: ${traj.size}</div>
            </div>
        `).join('');
    }

    updateCompareButton() {
        this.compareBtn.disabled = this.trajectories.length < 2;
    }

    async runComparison() {
        if (this.trajectories.length < 2) {
            this.showMessage('error', 'Please upload at least 2 trajectories');
            return;
        }

        this.compareBtn.disabled = true;
        this.showLoading();

        const requestData = {
            trajectory_ids: this.trajectories.map(t => t.id),
            analysis_type: this.analysisType.value,
            selection: this.selection.value,
            show_stats: this.showStats.checked,
            align: this.alignTrajectories.checked
        };

        try {
            const response = await fetch('/api/comparison/compare', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(requestData)
            });

            if (!response.ok) {
                throw new Error(`Comparison failed: ${response.statusText}`);
            }

            const result = await response.json();
            this.displayResults(result);
        } catch (error) {
            this.showMessage('error', `Comparison failed: ${error.message}`);
            this.resultsArea.innerHTML = `
                <div class="error-message">
                    <strong>Error:</strong> ${error.message}
                </div>
            `;
        } finally {
            this.compareBtn.disabled = false;
        }
    }

    showLoading() {
        this.resultsArea.innerHTML = `
            <div class="loading">
                <div class="spinner"></div>
            </div>
        `;
    }

    displayResults(result) {
        let html = '';

        // Display plots
        if (result.plots && result.plots.length > 0) {
            result.plots.forEach(plot => {
                html += `
                    <div class="plot-container">
                        <div class="plot-title">${plot.title}</div>
                        <img src="${plot.url}" alt="${plot.title}">
                    </div>
                `;
            });
        }

        // Display statistics
        if (result.statistics) {
            html += `
                <div class="plot-container">
                    <div class="plot-title">📊 Statistical Analysis</div>
                    ${this.renderStatistics(result.statistics)}
                </div>
            `;
        }

        // Display summary
        if (result.summary) {
            html += `
                <div class="plot-container">
                    <div class="plot-title">📝 Summary</div>
                    <p>${result.summary}</p>
                </div>
            `;
        }

        this.resultsArea.innerHTML = html;
    }

    renderStatistics(stats) {
        let html = '<table class="stats-table"><thead><tr>';

        // Table headers
        const headers = Object.keys(stats[Object.keys(stats)[0]]);
        html += '<th>Trajectory</th>';
        headers.forEach(header => {
            html += `<th>${header}</th>`;
        });
        html += '</tr></thead><tbody>';

        // Table rows
        Object.entries(stats).forEach(([trajName, values]) => {
            html += `<tr><td><strong>${trajName}</strong></td>`;
            headers.forEach(header => {
                const value = values[header];
                const formatted = typeof value === 'number' ? value.toFixed(3) : value;
                html += `<td>${formatted}</td>`;
            });
            html += '</tr>';
        });

        html += '</tbody></table>';
        return html;
    }

    showMessage(type, message) {
        const messageDiv = document.createElement('div');
        messageDiv.className = type === 'error' ? 'error-message' : 'success-message';
        messageDiv.textContent = message;

        this.resultsArea.insertBefore(messageDiv, this.resultsArea.firstChild);

        setTimeout(() => messageDiv.remove(), 5000);
    }

    formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';
        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
    }
}

// Initialize on page load
let comparisonManager;
document.addEventListener('DOMContentLoaded', () => {
    comparisonManager = new ComparisonManager();
});
