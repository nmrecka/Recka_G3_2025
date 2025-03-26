function listFilesInFolder() {

  var folderId = "";  // Replace with your folder's ID 
  var folder = DriveApp.getFolderById(folderId);
  var files = folder.getFiles();
  
  var data = [];
  // Add a header row if desired:
  data.push(["File Name", "Download URL"]);
  
  while (files.hasNext()) {
    var file = files.next();
    var fileId = file.getId();
    var downloadUrl = "https://lh3.googleusercontent.com/d/" + fileId;
    data.push([file.getName(), downloadUrl]);
  }
  
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var sheet = ss.getActiveSheet();
  sheet.clear(); // Clears existing content
  sheet.getRange(1, 1, data.length, data[0].length).setValues(data);
}